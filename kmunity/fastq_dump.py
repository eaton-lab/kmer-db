#!/usr/bin/env python

"""Call `fastq-dump {SRR} -O {tmpdir}` to download data.

"""

from typing import Optional
import tempfile
from subprocess import Popen, STDOUT, PIPE
from pathlib import Path
from loguru import logger
import requests
from utils import KmunityError


logger = logger.bind(name="kmunity")
SRA_VERSION = "3.0.5"


def get_sra_tools_linux(tmpdir: Path = None) -> Path:
    """Download sratools {version} to tmp dir.
    """
    tmpdir = Path(tmpdir if tmpdir is not None else tempfile.gettempdir())

    # url of recent tested version of SRA tools and destination
    url = (
        f"https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/{SRA_VERSION}"
        f"/sratoolkit.{SRA_VERSION}-ubuntu64.tar.gz"
    )
    tmptar = tmpdir / Path(url).name

    # download to tmp dir
    if not tmptar.exists():
        logger.info(f"Downloading sra toolkit {SRA_VERSION} to tmpdir")
        with open(tmptar, 'wb') as tz:
            res = requests.get(url, stream=True)
            res.raise_for_status()
            tz.write(res.raw.read())

    # extract archive into tmpdir
    sra_dir = tmpdir / f"sratoolkit.{SRA_VERSION}-ubuntu64/"
    if not sra_dir.exists():
        logger.debug("Extracting sra toolkit to tmpdir")
        cmd = ["tar", "xzvf", str(tmptar), "-C", str(tempfile.gettempdir())]
        with Popen(cmd, stderr=STDOUT, stdout=PIPE) as proc:
            comm = proc.communicate()
            if proc.returncode:
                raise KmunityError(comm[0].decode())
    return sra_dir


def get_fastq_dumped(
    sra_dir: Path,
    outpath: Path,
    srr: str,
    split: bool = False,
    subsample: Optional[int] = None,
) -> Path:
    """Return path to fastq file with reads for run SRR.

    """
    # build fastq command
    cmd = [
        str(sra_dir / "bin" / "fastq-dump"),
        str(srr),
        "-O", str(outpath),
    ]
    if split:
        cmd.extend(["--split-3"])
    if subsample:
        cmd.extend(["-X", str(int(subsample))])

    logger.info(f"Executing: {cmd}")
    with Popen(cmd, stderr=STDOUT, stdout=PIPE) as proc:
        comm = proc.communicate()
        if proc.returncode:
            raise KmunityError(comm[0].decode())
    return outpath / f"{srr}.fastq"


if __name__ == "__main__":
    sra_dir = get_sra_tools_linux()
    # srr = get_fastq_dumped(sra_dir, Path('/tmp'), 'SRR390728', 1000)
    srr = get_fastq_dumped(sra_dir, Path('/tmp'), 'SRR516269', split=0, subsample=10_000_000)
    print(srr)
