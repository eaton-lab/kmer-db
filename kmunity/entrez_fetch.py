#!/usr/bin/env python


"""Search NCBI API for SRR Illumina data not yet in database.
"""

from typing import List, Iterator, Dict, Any
# import time
import requests
import xml.etree.ElementTree as ET
from loguru import logger

logger = logger.bind(name="kmunity")
# import pandas as pd


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
        """Return ... for a single """
        xml = get_runinfo_batch_xml([srr])
        for rdict in iter_runinfo_from_xml(xml, min_bases):
            yield rdict
        # return next(iter_runinfo_from_xml(xml, min_bases))

    def iter_runinfo(self, database: str, min_bases: float = 1e9) -> Iterator[Dict[str, Any]]:
        """sample a random UID that is NOT YET IN the database.

        Fetches a _batch_ of entrez sra data, check for new sample not
        yet in database and which has sufficient data for kmer analysis.
        If no new samples in batch, then sample the next batch until
        a sample is found.
        """
        # Entrez search term for database WGS data
        term: str = MAMMALS_TERM if database == "mammals" else BIRDS_TERM

        # iterate over matches in the database
        for uids in self._iter_searched_uids(term):

            # skip any uids already in dataset
            # logger.info(uids)

            # fetch runinfo from NCBI as XML and parse it.
            xml = get_runinfo_batch_xml(uids)
            # logger.info(xml)

            # iterate over runs in chunk
            for run_dict in iter_runinfo_from_xml(xml, min_bases):
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


def iter_runinfo_from_xml(xml: str, min_bases: int) -> Iterator[Dict[str, Any]]:
    """Yield dicts w/ {srr, srs, organism, tax_id, bases} for each run
    in the xml.
    """
    tree = ET.fromstring(xml)
    for exppack in tree:
        for exp in exppack:
            if exp.tag == "RUN_SET":
                for run in exp:
                    srr_accession = run.attrib["accession"]
                    srr_total_bases = int(run.attrib["total_bases"])

                    if srr_accession in EXCLUDE_TAXIDS:
                        logger.debug("skipping {srr_accession}; already done.")
                        continue

                    if not srr_total_bases > min_bases:
                        logger.debug("skipping {srr_accession}; too few bases.")
                        continue

                    for pool in run:
                        if pool.tag == "Pool":
                            for member in pool:
                                srs_accession = member.attrib["accession"]
                                organism = member.attrib["organism"]
                                tax_id = int(member.attrib["tax_id"])

                                yield {
                                    "SRR": srr_accession,
                                    "SRS": srs_accession,
                                    "organism": organism,
                                    "tax_id": tax_id,
                                    "bases": srr_total_bases,
                                }


if __name__ == "__main__":

    e = EntrezFetch()
    # get a specific Run
    e.get_runinfo("SRR7811753")
    # get ordered matches in a database
    birds = e.iter_runinfo("birds")
    for i in range(10):
        print(next(birds))