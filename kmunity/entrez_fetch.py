#!/usr/bin/env python

"""Search NCBI API for SRR Illumina data not yet in database.

"""

from typing import List, Iterator, Dict, Any
from pathlib import Path
import xml.etree.ElementTree as ET
import requests
from loguru import logger
import pandas as pd

logger = logger.bind(name="kmunity")


MAMMALS_TERM = """
(((((((Mammalia[Organism])
AND "public"[Access])
AND "illumina"[Platform])
AND "wgs"[Strategy])
AND "genomic"[Source])
AND "filetype fastq"[Filter])
AND "sra nuccore wgs"[Filter])
AND "strategy whole genome sequencing"[Filter]
""".replace("\n", "")


BIRDS_TERM = """
(((((((Aves[Organism])
AND "public"[Access])
AND "illumina"[Platform])
AND "wgs"[Strategy])
AND "genomic"[Source])
AND "filetype fastq"[Filter])
AND "sra nuccore wgs"[Filter])
AND "strategy whole genome sequencing"[Filter]
""".replace("\n", "")


# exclude Homo sapiens, Mus musculus, Bos_x_
EXCLUDE_TAXIDS = [9606, 10090, 9615]


class EntrezFetch:
    """Search NCBI for a match to a search term from a supported database.

    Example
    -------
    >>> e = EntrezFetch()

    >>> # get a specific Run
    >>> e.get_runinfo("SRR7811753")

    >>> # get ordered matches in a database
    >>> birds = e.iter_runinfo("birds")
    >>> for i in range(10):
    >>>     print(next(birds))
    """
    def __init__(self, database: str):
        # Entrez search term for database WGS data
        assert database in ("mammals", "birds")
        self.database = database
        self.term = MAMMALS_TERM if database == "mammals" else BIRDS_TERM
        self.data = pd.read_csv(Path(f"../{database}") / "database.csv", index_col=0)

    def _iter_searched_uids(self, term, chunksize: int = 20) -> Iterator[List[int]]:
        """Yield a chunk of uids matching the search term.
        """
        # get number of NCBI Entrez records matching search term
        nuids = count_uids(term)

        # yield [chunksize] uids at a time until exhausted
        for chunk in range(0, nuids, chunksize):
            uids = get_uids(term, retstart=chunk, retmax=chunk + chunksize)
            yield uids

    def get_runinfo(self, srr: str = None, min_bases: float = 1e9) -> Dict[str, Any]:
        """Return runinfo for a single SRR if it passes min-bases filter"""
        xml = get_runinfo_batch_xml([srr])
        return list(iter_runinfo_from_xml(xml, min_bases))

    def iter_runinfo(self, min_bases: float = 1e9) -> Iterator[Dict[str, Any]]:
        """Generator of RunInfo dicts for SRR runs not yet in database.

        Fetches a _batch_ of entrez sra data, check for new sample not
        yet in database and which has sufficient data for kmer analysis.
        If no new samples in batch, then sample the next batch until
        a sample is found.
        """
        # iterate over matches in the database
        for uids in self._iter_searched_uids(self.term):

            # fetch runinfo from NCBI as XML and parse it.
            xml = get_runinfo_batch_xml(uids)

            # iterate over runs in chunk
            for run_dict in iter_filtered_runinfo(xml, min_bases, self.data):
                yield run_dict


def count_uids(term: str) -> int:
    """Search NCBI nucleotide database and return N sequence matches"""

    # make a request to esearch
    res = requests.get(
        url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
        params={
            "db": "sra",
            "term": term,
            "retmode": "text",
            "tool": "kmunity",
            "email": "research@univerity.edu"
        },
    )

    # parse the xml output
    count = res.text.split("<Count>")[1].split("</Count>")[0].strip()
    return int(count)


def get_uids(term: str, retstart: int = 0, retmax: int = 200) -> List[int]:
    """Search NCBI nucleotide database and return UIDS for N matches to
    search term. The max searches allowed SEEMS to be around 5000.
    """
    # make a request to esearch
    res = requests.get(
        url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
        params={
            "db": "sra",
            "term": term,
            "retmode": "text",
            "retstart": str(retstart),
            "retmax": str(retmax),
            "tool": "kmunity",
            "email": "research@univerity.edu"
        },
    )

    # parse the xml output
    count = res.text.split("<Count>")[1].split("</Count>")[0].strip()

    # if nothing found then bail out
    if not int(count):
        raise ValueError("No UIDs found")

    # return the list of UIDs from xml, e.g.: <Id>27123219</Id>
    uids = []
    ids = res.text.split("<IdList>")[1].split("</IdList>")[0].strip()
    for item in ids.split("\n"):
        uids.append(item[4:-5])
    return uids


def get_runinfo_batch_xml(uids: List[int]) -> str:
    """Search NCBI for metadata (runinfo) by supplying UIDs. The max
    allowed searches for this seems to be small, like 20 or so.
    """
    res = requests.get(
        url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
        params={
            "db": "sra",
            "id": ",".join(uids),
            "retmode": "text",
            "tool": "kmunity",
            "email": "research@univerity.edu"
        }
    )
    return res.text


def iter_runinfo_from_xml(xml: str) -> Iterator[Dict[str, Any]]:
    """Yield dicts of {srr, srs, organism, tax_id, bases} for each run in xml.
    """
    # iterate over RUN items in XML to extract sample data
    tree = ET.fromstring(xml)
    for run in tree.findall(".//RUN"):

        # get RUN (SRR) accession
        srr_accession = run.attrib["accession"]

        # get BASES
        srr_total_bases = int(run.attrib["total_bases"])

        # iterate over SAMPLES (SRS) in this RUN
        for pool in run.findall(".//Pool"):

            samn_accession = pool.find('.//EXTERNAL_ID').text
            for member in pool:

                # Sample data
                srs_accession = member.attrib["accession"]
                organism = member.attrib["organism"]
                tax_id = int(member.attrib["tax_id"])

                # runinfo dict
                yield {
                    "SRR": srr_accession,
                    "SRS": srs_accession,
                    "SAMN": samn_accession,
                    "organism": organism,
                    "tax_id": tax_id,
                    "bases": srr_total_bases,
                }


def iter_filtered_runinfo(
    xml: str,
    min_bases: int,
    data: pd.DataFrame,
) -> Iterator[Dict[str, Any]]:
    """Return runinfo dicts passing filters

    (1): Sufficient number of bases (min_bases)
    (2): Not already in database (data, filter)

    For example, only return if tax_id, SRR, SRS, or SAMN is not already
    in database.
    """
    for rdict in iter_runinfo_from_xml(xml):

        if rdict["bases"] < min_bases:
            logger.info(f"skipping {rdict['SRR']}; too few bases")
            continue

        if rdict["SRR"] in data.SRR:
            logger.info(f"skipping {rdict['SRR']}; already in database")
            continue

        if rdict["SRS"] in data.SRS:
            logger.info(f"skipping {rdict['SRA']}; SRS ID ({rdict['SRS']}) already in database")
            continue

        if rdict["SAMN"] in data.SAMN:
            logger.info(f"skipping {rdict['SRR']}; SAMN ID ({rdict['SAMN']}) already in database")
            continue

        if rdict["tax_id"] in data.tax_id:
            logger.info(f"skipping {rdict['SRR']}; TAX ID ({rdict['tax_id']}) already in database")
            continue

        # enter to dataframe and yield as dict
        data.loc[rdict["SRS"]] = [rdict.get(i) for i in data.columns]
        yield rdict


if __name__ == "__main__":

    import kmunity
    kmunity.set_log_level("DEBUG")

    e = EntrezFetch("birds")

    # get a specific Run
    # r = e.get_runinfo("SRR7811753")
    # print(r)

    # get ordered matches in a database
    birds = e.iter_runinfo(min_bases=2e9)
    for i in range(20):
        print(next(birds))
    print(e.data)