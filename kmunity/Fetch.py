#!/usr/bin/env python


"""
Search NCBI API for SRR Illumina data not yet in database.
"""

import time
import requests
import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np


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



class SearchSRR(object):
    """
    Search NCBI for a specific SRR match.
    """
    def __init__(self, srr):
        self.srr = srr


    def run(self, sleep=10):
        # get internal ids for this search. A list of matches.
        self.uids = get_uids(self.srr)
        time.sleep(sleep)

        # get runinfo for all files associated with this run
        self.xml = get_runinfo(self.uids[:10])
        time.sleep(sleep)

        # extract srr result for the .sra release
        self.res = parse_runinfo(self.xml)

        # if the extraction worked, then continue
        self.bio, self.org, self.tax, self.bases, self.run = self.res[0]



class SearchTerm(object):
    """
    Search NCBI for a match to a search term from a supported database.
    """
    def __init__(self, db):
        self.db = db
        self.term = MAMMALS_TERM  # if self.db = "mammals" else "birds")


    def run(self, sleep=5.5):
        """
        sample a random UID that is NOT YET IN the database.
        """
        self.nuids = count_uids(self.term)
        self.uids = get_uids(self.term, retstart=0, retmax=self.nuids)

        # continue until we find a match.
        self.rinfo = None
        ntries = 0
        while not self.rinfo:
            self.uid = [self.uids[np.random.randint(self.nuids)]]
            self.xml = get_runinfo(self.uid)
            self.rinfo = parse_runinfo(self.xml, mincov_gb=0)
            ntries += 1
            time.sleep(sleep)

        # if the extraction worked, then continue
        # print(self.uid, ntries, self.rinfo)
        self.bio, self.org, self.tax, self.bases, self.run = self.rinfo[0]



def count_uids(term):
    "Search NCBI nucleotide database and return N sequence matches"

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




def get_uids(term, retstart=0, retmax=200):
    """
    Search NCBI nucleotide database and return UIDS for N matches to search term. 
    The max searches allowed SEEMS to be around 5000.
    """

    # make a request to esearch 
    res = requests.get(
        url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi", 
        params={
            "db": "sra",
            "term": term,
            "retmode": "text",
            "retstart": retstart,
            "retmax": retmax,
            "tool": "kmunity", 
            "email": "research@univerity.edu"
            },
        )

    # parse the xml output
    count = res.text.split("<Count>")[1].split("</Count>")[0].strip()

    # if nothing found then bail out
    if not int(count):
        raise ValueError("No UIDs found")

    # return the list of UIDs
    uids = []
    ids = res.text.split("<IdList>")[1].split("</IdList>")[0].strip()
    for item in ids.split("\n\t"):
        uids.append(item[4:-5])
    return uids    




def get_runinfo(uids):
    """
    Search NCBI for metadata (runinfo) by supplying UIDs. The max allowed 
    searches for this seems to be small, like 20 or so.
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




def parse_runinfo(xml, mincov_gb=0, exclude_taxids=EXCLUDE_TAXIDS):
    """
    Sift through UID matches to find hits with WGS data of sufficient amounts.
    """
    tree = ET.fromstring(xml)
    data = []
    tax_id = None
    for exppack in tree:
        for exp in exppack:
            if exp.tag == "RUN_SET":
                for run in exp:
                    for pool in run:
                        if pool.tag == "Pool":
                            for member in pool:
                                # print(member.tag, member.attrib)
                                accession = member.attrib["accession"]
                                organism = member.attrib["organism"]
                                tax_id = int(member.attrib["tax_id"])
                                gbases = int(int(member.attrib["bases"]) / 1e9)

                        if pool.tag == "SRAFiles":
                            for member in pool:
                                if member.attrib["semantic_name"] == "run":
                                    # print(member.tag,  member.attrib)
                                    SRR = member.attrib["filename"]

        # exclude some from table
        if tax_id:
            if gbases > mincov_gb:
                if tax_id not in exclude_taxids:
                    data.append([accession, organism, tax_id, gbases, SRR])
                else:
                    pass
                    # print("dropped: organism in EXCLUDED_TAXON set.")
            else:
                pass
                # print("dropped: too little data.")
    return data



if __name__ == "__main__":
    pass