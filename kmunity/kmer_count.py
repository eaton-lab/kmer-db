#!/usr/bin/env python

"""Count kmers from fastq data with kmc

"""

from typing import Dict
import tempfile
from subprocess import Popen, STDOUT, PIPE
from pathlib import Path
from loguru import logger
import requests
from utils import KmunityError


logger = logger.bind(name="kmunity")


def get_gce_tools_linux(tmpdir: Path = None) -> Dict[str, Path]:
    """Return dict mapping Path to gce and kmerfreq binaries.
    """
    tmpdir = Path(tmpdir if tmpdir is not None else tempfile.gettempdir())

    # url of recent tested version of SRA tools and destination
    urls = [
        "https://github.com/fanagislab/GCE/raw/master/gce-1.0.2/gce",
        "https://github.com/fanagislab/GCE/raw/master/gce-1.0.2/kmerfreq",
    ]

    binaries = {}
    for url in urls:
        res = requests.get(url, allow_redirects=True)
        exe = Path(url).name
        outbin = tmpdir / exe

        if not outbin.exists():
            logger.info(f"Downloading {exe} binary from {url}")
            with open(outbin, 'wb') as out:
                out.write(res.content)

                # ensure it is executable
                cmd = ['chmod', '+x', outbin]
                with Popen(cmd, stderr=STDOUT, stdout=PIPE) as proc:
                    out = proc.communicate()
        binaries[exe] = outbin
    return binaries


def get_kmers_from_fastq(
    kmerfreq: Path,
    fastq: Path,
    outpath: Path,
    kmer_size: int = 17,
    threads: int = 4,
) -> Path:
    """

    """
    # call the tool
    cmd = [
        str(kmerfreq),
        "-k", str(kmer_size),
        "-t", str(threads),
        "-p", str(outpath),
        str(fastq)
    ]

    # do not log local files paths
    logger.info("Executing: {kmerfreq} ... {fastq}")
    logger.debug(f"Executing: {' '.join(cmd)}")
    with Popen(cmd, stderr=STDOUT, stdout=PIPE) as proc:
        out = proc.communicate()
        if proc.returncode:
            raise KmunityError(out[0].decode())

        # log head results file
        # resfile = f"{fastq}.srr" + ".kmer.freq.stat"
        # logger.success("Kmer counts complete: {}".format(resfile))


if __name__ == "__main__":

    tools = get_gce_tools_linux()
    print(tools)
