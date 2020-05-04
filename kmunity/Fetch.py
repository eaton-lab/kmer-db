#!/usr/bin/env python


"""
Search NCBI API for SRR Illumina data not yet in database.
"""

import requests
import xml.etree.ElementTree as ET
import pandas as pd
import time



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



class Search_SRR:
    def __init__(self, srr):
        self.srr = srr


    def run(self):
        # get internal ids for this search
        self.uids = get_uids(self.srr)
        time.sleep(10)

        # get runinfo for all files associated with this run
        self.xml = get_runinfo(self.uids)
        time.sleep(10)        

        # extract srr result for the .sra release
        self.res = parse_runinfo(self.xml)
        self.bio, self.org, self.tax, self.bases, self.run = self.res[0]



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




def get_uids(term, retstart=0, retmax=20):
    "Search NCBI nucleotide database and return N sequence matches"

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
    for item in ids.split("\n"):
        uids.append(item[4:-5])
    return uids    



def get_runinfo(uids):
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




def parse_runinfo(xml, mincov_gb=5, exclude_taxids=EXCLUDE_TAXIDS):
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
                                bases = int(int(member.attrib["bases"]) / 1e9)

                        if pool.tag == "SRAFiles":
                            for member in pool:
                                if member.attrib["semantic_name"] == "run":
                                    # print(member.tag,  member.attrib)
                                    SRR = member.attrib["filename"]

        # exclude some from table
        if tax_id:
            if bases > 5:
                if tax_id not in exclude_taxids:
                    data.append([accession, organism, tax_id, bases, SRR])
    return data
