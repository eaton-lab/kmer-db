#!/usr/bin/env python

"""
Download fastq files to calculate genome size and heterozygosity with gce
"""

from __future__ import print_function

import os
import sys
import glob
import tempfile
import subprocess as sps

import requests
import pandas as pd
import kmunity
from loguru import logger

from .Fetch import Search_SRR



class Kmunity:
    """
    By default it looks for software in your PATH. 

    Parameters
    ----------
    srr: (str or None)
        NCBI Run ID, e.g., SRR7811753. If None then a sample that does not yet
        exist in the database will be fetched from NCBI for the selected 
        database (db).

    db: (str)
        The database to which a sample belongs (currently supported includes
        mammals, birds, plants). If an 'srr' is entered it will be checked
        for appropriate inclusion in the database.

    workdir: (str):
        A temporary directory in which the downloaded data can be stored. This
        should have plenty of space available. kmunity will remove all files 
        downloaded to this location after finishing. 

    repo: (str)
        The kmunity git repository where current results are stored and where
        new results will be organized. Default is "./kmunity".
    """
    def __init__(self, srr=None, db="mammals", workdir="/tmp", repo="./kmunity", **kwargs):

        # store args
        self.srr = srr  # (srr if srr else "SRRXXXYYY")
        self.db = db
        self.data = None

        # expand i/o file paths
        self.repo = os.path.realpath(os.path.expanduser(repo))
        self.csv = os.path.join(self.repo, self.db, "database.csv")
        self.workdir = os.path.realpath(os.path.expanduser(workdir))
        self.srrdir = os.path.join(self.workdir, self.srr)
        self.logdir = os.path.join(self.repo, self.db, "logfiles") 
        self.logfile = os.path.join(self.logdir, "{}.log".format(self.srr))
        self._logger_set()

        # check kwargs: e.g., user-supplied binary paths
        self.binaries = {}

        # allow kwargs to overwrite binary paths
        for key in kwargs:
            if key in self.binaries:
                self.binaries[key] = kwargs[key]

        # config setup only
        if kwargs.get("config"):
            self._get_binary()
            self._vdb_config()

        # run checks on existing results, paths and binaries.
        else:
            self._get_usergh()
            self._get_binary()
            self._path_check()
            self._query_ncbi()



    def _logger_set(self):
        """
        Configure Loguru to log to stdout and logfile.
        """
        # add stdout logger
        config = {
            "handlers": [
                {
                    "sink": sys.stdout, 
                    "format": (
                        "{time:YYYY-MM-DD-hh:mm} | "
                        "<cyan>{function}</cyan> | "
                        "<level>{message}</level>"
                    ),
                    "level": "DEBUG",
                    },
                {
                    "sink": self.logfile,                   
                    "format": "{time:YYYY-MM-DD} | {function} | {message}",
                    "level": "INFO",
                    }
            ],
            "extra": {"user": "deren"}
        }
        logger.configure(**config)
        logger.enable("kmunity")

        # if logfile exists then report it
        if os.path.exists(self.logfile):

            if 0:
                # if completed file found
                logger.warning(
                    "A completed logfile exists for this accession from {}.\n"
                    "Look in the database file for results."
                )

            else:
                # if logfile is unfinished (not completed run) then remove file
                logger.debug("")
                logger.debug("LOG FILE")
                logger.debug(
                    'Clearing previous unfinished run of {}.'
                    .format(self.srr))
                open(self.logfile, 'w').close()
                logger.debug("")



    def _path_check(self):
        """
        Store path locations and check for existing results.
        """
        # ensure repo path is correct
        assert os.path.exists(self.repo), (
            "'repo' path not found: {}".format(self.repo))
        assert os.path.exists(self.csv), (
            "'repo' path does not point to local kmunity repo: {}"
            .format(self.csv))

        # ensure workdir and logdir exist
        for dirname in [self.workdir, self.logdir, self.srrdir]:
            if not os.path.exists(dirname):
                os.makedirs(dirname)

        # load existing database
        self.data = pd.read_csv(self.csv)
        logger.debug("LOCAL PATHS")
        logger.debug("workdir: {}".format(self.workdir))
        logger.debug("srrdir: {}".format(self.srrdir))        
        logger.debug("logfile: {}".format(self.logfile))
        logger.debug("database: {}".format(self.csv))
        logger.debug("")        


        # curdb = pd.read_csv(self.csv)
        # if not os.path.exists(self.kwargs['outdir']):
            # os.makedirs(self.kwargs['outdir'])



    def _query_ncbi(self):
        """
        Query NCBI to get info for user selected SRR accession, OR, to find
        an SRR accession from the selected database that is not yet in the 
        kmunity database.
        """
        # No user SRR supplied
        if self.srr is None:
            logger.error("No SRR option not yet supported.")
            raise NotImplementedError("No SRR option not yet supported.")

        # user SRR supplied
        else:
            # print a log header with info
            logger.warning("NCBI QUERY")
            self.query = Search_SRR(self.srr)
            self.query.run()
            logger.info("database: {}".format(self.db))
            logger.info("organism: {}".format(self.query.org))
            logger.info("taxid: {}".format(self.query.tax))
            logger.info("biosample: {}".format(self.query.bio))
            logger.info("run: {}".format(self.query.run))
            logger.info("size (Gb): {}".format(self.query.bases))
            # logger.info("status: {}".format(bool(self.srr in self.data.Run))
            logger.info("")



    def _get_usergh(self):
        # get user
        user = "anonymous"
        try:
            cmd = ["git", "config", "--get", "user.name"]
            proc = sps.Popen(cmd, stderr=sps.STDOUT, stdout=sps.PIPE)
            user = proc.communicate()[0].decode().strip()
        except Exception:
            pass
        logger.warning("CONTRIBUTOR")
        logger.info("GitHub user: {}".format((user if user else "unknown")))
        if not user:
            logger.debug(
                "tip: set username with: git config --global user.name 'Hotdog'")
        logger.info("")



    def _get_binary(self):
        """
        Always pulls in fixed versions of gce and kmerfreq to srrdir.
        """
        bin_gce = os.path.join(tempfile.gettempdir(), "gce")
        bin_kme = os.path.join(tempfile.gettempdir(), "kmerfreq")
        bin_pre = os.path.join(
            tempfile.gettempdir(), 
            "sratoolkit.2.10.8-ubuntu64/"
            "bin", "prefetch")
        bin_fas = os.path.join(
            tempfile.gettempdir(), 
            "sratoolkit.2.10.8-ubuntu64/"
            "bin", "fasterq-dump")
        bin_vdb = os.path.join(
            tempfile.gettempdir(), 
            "sratoolkit.2.10.8-ubuntu64/"
            "bin", "vdb-config")

        if not os.path.exists(bin_gce):
            logger.debug("local gce not found.")
            self._dl_gce_tmp()
            if not os.path.exists(bin_gce):
                logger.error("gce tools download failed.")
        self.binaries["gce"] = bin_gce
        self.binaries["kmerfreq"] = bin_kme

        if not os.path.exists(bin_pre):
            logger.debug("local sratools >=2.10.8 not found.")            
            self._dl_sra_tmp()
            if not os.path.exists(bin_pre):
                logger.error("sratools download failed.")
        self.binaries["prefetch"] = bin_pre
        self.binaries["fasterq-dump"] = bin_fas
        self.binaries["vdb-config"] = bin_vdb

        # print software versions
        logger.warning("VERSIONS")
        logger.info("kmunity: {}".format(kmunity.__version__))
        logger.info("sra-tools: {}".format(self._x_prefetch(True)))
        logger.info("kmerfreq: {}".format(self._x_kmerfreq(True)))
        logger.info("gce: {}".format(self._x_call_gce(True)))
        logger.info("")



    def _dl_gce_tmp(self):
        """
        Downloads gce and kmerfreq binaries, checks +x, and locates to tmp.
        """
        # pull gce & kmer executables to workdir
        logger.debug("Downloading gce and kmerfreq to /tmp")
        gce_url = "https://github.com/fanagislab/GCE/raw/master/gce-1.0.2/gce"
        kme_url = "https://github.com/fanagislab/GCE/raw/master/gce-1.0.2/kmerfreq"
        for url in [gce_url, kme_url]:
            res = requests.get(url, allow_redirects=True)
            exe = os.path.basename(url)
            outbin = os.path.join(tempfile.gettempdir(), exe)
            with open(outbin, 'wb') as out:
                out.write(res.content)
                self.binaries[exe] = outbin

            # ensure it is executable
            cmd = ['chmod', '+x', outbin]
            proc = sps.Popen(cmd, stderr=sps.STDOUT, stdout=sps.PIPE)
            out = proc.communicate()



    def _dl_sra_tmp(self):
        """
        Downloads sratools (linux), checks +x, and locates to tmp.
        """
        logger.debug("Downloading sratoolkit.2.10.8 to /tmp")
        # pull in version 2.10.5 of sra-tools (ONLY THIS VERSION WORKS)
        # https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.8/
        url = "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.8/sratoolkit.2.10.8-ubuntu64.tar.gz"
        tmptar = os.path.join(tempfile.gettempdir(), os.path.basename(url))
        res = requests.get(url, stream=True)
        res.raise_for_status()
        with open(tmptar, 'wb') as tz:
            # for chunk in res.iter_content(chunk_size=8192):
            tz.write(res.raw.read()) 

        # decompress tar file 
        logger.debug("Extracting sratoolkit.2.10.8 in /tmp")
        logger.debug("")
        cmd = ["tar", "xzvf", tmptar, "-C", tempfile.gettempdir()]
        proc = sps.Popen(cmd, stderr=sps.STDOUT, stdout=sps.PIPE)
        comm = proc.communicate()
        if proc.returncode:
            print(comm[0].decode())



    def _x_prefetch(self, version_only=False):
        """
        - if version: Check for prefetch tool from sratools by calling -V
        - Calls prefetch on self.srr and reports .sra download success.
        """
        if version_only:
            # print the version
            proc = sps.Popen(
                [self.binaries["prefetch"], "-V"], 
                stderr=sps.STDOUT, 
                stdout=sps.PIPE,
            )
            out = proc.communicate()
            if proc.returncode:
                logger.error("prefetch tool not found.")
            return out[0].decode().split()[-1]

        # log the command used to prefetch
        cmd = [
            self.binaries["prefetch"], self.srr, 
            "-O", self.srrdir,
            "-X", str(int(1e9)),
        ]
        logger.info("Executing: {prefetch} {srr} -O {srrdir} -X 1000000000")
        logger.debug("Executing: {}".format(" ".join(cmd)))

        # call execute        
        proc = sps.Popen(cmd, stderr=sps.STDOUT, stdout=sps.PIPE)
        out = proc.communicate()
        if proc.returncode:
            logger.error("Error: {}".format(out[0].decode()))
            raise Exception(out[0].decode())

        # show file size result
        srafile = os.path.join(self.srrdir, self.srr + ".sra")
        if not os.path.exists(srafile):
            logger.error("Prefetch failed, no sra file found.")
            raise IOError("Prefetch failed, no sra file found.")
        size = os.path.getsize(srafile)
        size = round(size / 1e9, 2)
        logger.success("Downloaded {} ({} Gb)".format(self.srr + ".sra", size))



    def _x_fasterqd(self, version_only=False):
        """
        
        """
        if version_only:
            # print the version
            proc = sps.Popen(
                [self.binaries["fasterq-dump"], "-V"],
                stderr=sps.STDOUT, 
                stdout=sps.PIPE,
            )
            out = proc.communicate()[0].decode().split()[-1]
            return out

        # commands to the logger
        cmd = [
            self.binaries["fasterq-dump"], self.srr, 
            "-O", self.srrdir,
            "-t", self.srrdir,
        ]
        null = "{fasterq-dump} {srr} -O {workdir}/{srr} -t {workdir}/{srr}"
        logger.info("Executing: {}".format(null))
        logger.debug("Executing: {}".format(" ".join(cmd)))

        # call the tool
        proc = sps.Popen(cmd, stderr=sps.STDOUT, stdout=sps.PIPE)
        out = proc.communicate()
        if proc.returncode:
            logger.error("Failed: {}".format(out[0].decode()))
            logger.error(
                "If you encountered a disk full error but believe you have "
                "sufficient space available in your working dir then the "
                "is being caused by the insane behavior of the sra-tools "
                "package which hides large tmp file in obscure places. "
                "You can turn off this behavior by running 'vdb-config -i'."
                "Turn off the 'enable local file caching' option."
                )
            raise TypeError(out[0].decode())

        # write a tmp SRR.lib file
        libfile = os.path.join(self.srrdir, "{}_files.lib".format(self.srr))
        fastqs = glob.glob(os.path.join(self.srrdir, "*.fastq"))
        with open(libfile, 'w') as out:
            out.write("\n".join(fastqs))

        # show file size result
        f1 = self.srr + "_1.fastq"
        fastq1 = os.path.join(self.srrdir, f1)
        size1 = os.path.getsize(fastq1)
        size1 = round(size1 / 1e9, 2)

        f2 = self.srr + "_2.fastq"
        fastq2 = os.path.join(self.srrdir, f2)
        size2 = os.path.getsize(fastq2)
        size2 = round(size2 / 1e9, 2)

        logger.success("Fastq dumped {} ({} Gb)".format(f1, size1))
        logger.success("Fastq dumped {} ({} Gb)".format(f2, size2))



    def _x_kmerfreq(self, version_only=False):

        if version_only:
            proc = sps.Popen(
                [self.binaries["kmerfreq"], "-h"],
                stderr=sps.STDOUT, 
                stdout=sps.PIPE,
            )
            out = proc.communicate()[0].decode().split("\n")
            out = [i for i in out if "Version" in i][0]
            vers = out.split()[-1]
            return vers


        # call the tool
        lib = os.path.join(self.srrdir, "{}_files.lib".format(self.srr))
        cmd = [
            self.binaries["kmerfreq"],
            "-k", "17",
            "-t", "4",
            "-p", os.path.join(self.srrdir, self.srr),
            lib,
        ]

        # do not log local files paths
        null = "{kmerfreq} -k 17 -t 4 -p {srrdir}/{srr} {srrdir}/{srr}_files.lib"
        logger.info("Executing: {}".format(null))
        logger.debug("Executing: {}".format(" ".join(cmd)))
        proc = sps.Popen(cmd, stderr=sps.STDOUT, stdout=sps.PIPE)
        out = proc.communicate()
        if proc.returncode:
            raise Exception(out[0].decode())

        # log head results file
        resfile = self.srr + ".kmer.freq.stat"
        logger.success("Kmer counts complete: {}".format(resfile))

        # with open(os.path.join(self.srrdir, resfile), 'r') as indata:
            # logger.info("FILE CONTENTS:\n" + "".join(indata.read()))



    def _vdb_config(self):
        """
        Check whether vdb has been config'd. If not, warn user and exit.
         Check whether vdb-config has cache true, and exit if yes.
        """
        # user config file is F&$&*$# FORCED to be in home by sra-tools.
        logger.warning("CONFIG")
        logger.info(
            "To finish config sra-tools *requires* you run vdb-config once.")
        logger.info(
            "I highly recommend turning off 'enable local file caching'.")
        logger.info(
            "{} -i"
            .format(self.binaries["vdb-config"]))
        logger.debug(
            "This will create a sra-tools config file in {}"
            .format(os.path.expanduser("~/.ncbi/")))

        logger.debug("")
        # user config file is F&$&*$# FORCED to be in home by sra-tools.




    # def _xcutadapt(self, version_only=False):

    #     if version_only:
    #         # print the version
    #         proc = sps.Popen(
    #             [self.binaries["cutadapt"], "-v"],
    #             stderr=sps.STDOUT, 
    #             stdout=sps.PIPE,
    #         )
    #         out = proc.communicate()[0].decode().split("\n")
    #         out = [i for i in out if "Version" in i][0]
    #         vers = out.split()[-1]
    #         return vers

    #     # build new command
    #     cmd = [
    #         self.binaries["cutadapt"], 
    #         "-j", "4",
    #         "-o", self.fastq1, 
    #         "-p", self.fastq2,
    #         self.clean1, 
    #         self.clean2,
    #     ]





    def _x_call_gce(self, version_only=False):

        if version_only:
            # print the version
            proc = sps.Popen(
                [self.binaries["gce"], "-V"],
                stderr=sps.STDOUT, 
                stdout=sps.PIPE,
            )
            out = proc.communicate()[0].decode().split("\n")
            out = [i for i in out if "Version" in i][0]
            vers = out.split()[-1]
            return vers

        # prerun commands 
        resfile = os.path.join(self.srrdir, self.srr + ".kmer.freq.stat")
        cmd1 = ['cat', resfile]
        cmd2 = ['grep', '#Kmer indivdual number']  # (sic)
        null = "cat {srrdir}/{srr}.kmer.freq.stat | grep '#Kmer indiv'"
        logger.info("Executing {}:".format(null))
        logger.debug("Executing: {}".format(" ".join(cmd1)))
        proc1 = sps.Popen(cmd1, stderr=sps.STDOUT, stdout=sps.PIPE)
        proc2 = sps.Popen(cmd2, stdin=proc1.stdout, stdout=sps.PIPE)
        out = proc2.communicate()
        if proc2.returncode:
            logger.error(out[0].decode())

        # store ikmer result
        ikmer = out[0].decode().strip().split()[-1]
        logger.success("Kmer individual number: {}".format(ikmer))

        # write kmer 2 column sub result file
        res2col = resfile + ".2colum"
        logger.info("Parsing 2-columns file to: {}".format(res2col))
        arr = pd.read_csv(resfile, skiprows=7, sep="\t", header=None)
        arr = arr.iloc[:, :2]
        arr.to_csv(res2col, index=False, sep="\t", header=None)

        # Run in homozygous mode
        logger.info("Running 'gce' in homozygous mode to estimate coverage")
        null = "{gce} -g " + ikmer + " -f {res.2col}"
        logger.info("Executing: {}".format(null))
        cmd = [
            self.binaries["gce"],
            "-g", ikmer,
            "-f", res2col,
        ]
        logger.debug(" ".join(cmd))
        proc = sps.Popen(cmd, stderr=sps.STDOUT, stdout=sps.PIPE)
        self.gce1out = proc.communicate()
        if proc.returncode:
            logger.error(self.gce1out[0].decode())

        # write to a tmp file
        parse = self.gce1out[0].decode().split("Final estimation table:")
        parse = parse[-1].strip().split("\n")
        headers, data = parse[:2]
        headers = headers.strip().split("\t")
        data = data.strip().split("\t")
        self.h0dict = {}
        for i, j in zip(headers, data):
            self.h0dict[i] = j
            logger.success("GCE H0 {}: {}".format(j))

        # Run in heterozygous mode
        logger.info("Running 'gce' in heterozygous mode.")
        null = "{gce} -g " + ikmer + " -f {res.2col} -H 1 -c {coverage}"
        logger.info("Executing: {}".format(null))
        cmd = [
            self.binaries["gce"],
            "-g", ikmer,
            "-f", res2col,
            "-H", "1", 
            "-c", str(int(float(round(self.h0dict["coverage_depth"]))))
        ]
        logger.debug(" ".join(cmd))
        proc = sps.Popen(cmd, stderr=sps.STDOUT, stdout=sps.PIPE)
        self.gce2out = proc.communicate()
        if proc.returncode:
            logger.error(self.gce2out[0].decode())

        # write to a tmp file
        parse = self.gce2out[0].decode().split("Final estimation table:")
        parse = coverage[-1].strip().split("\n")
        logger.success("GCE genome size: {}".format(coverage))
        logger.success("GCE heterozygosity: {}".format(coverage))
        logger.success("GCE coverage depth: {}".format(coverage))




    def parse_results(self):
        """
        Copy results from logfile to the CSV database table.
        [uid, organism, taxid, biosample, run, data-size, k, g-size, cov, het, filter-opts]
        """
        pass


    def _clean_work(self):
        """
        Remove any temp files.
        """
        logger.info("Removing temp files in {workdir}")
        logger.debug("Removing: {}".format(self.srrdir))




    # THIS DOESN'T WORK, VDB-CONFIG CANNOT DISABLE CACHE NON-INTERACTIVELY
    # AS FAR AS I CAN TELL. UGH. NEED TO ASK USERS TO DO IT WITH -I.
    # def _set_vdbcfg(self, reset=False):
    #     # reset to initial value
    #     if reset:
    #         cmd0 = [
    #             "vdb-config", "--set", "cache-disabled:{}"
    #             .format(self.init_cache)]
    #         out = sps.Popen(cmd0).communicate()
    #         logger.debug(
    #             "{vdb-config} --set cache-disabled:{}"
    #             .format(self.init_cache))

    #     # get and store the initial cache setting
    #     cmd0 = [self.binaries["vdb-config"], "-a"]
    #     cmd1 = ["grep", "cache"]
    #     proc0 = sps.Popen(cmd0, stderr=sps.STDOUT, stdout=sps.PIPE)
    #     proc1 = sps.Popen(cmd1, stdin=proc0.stdout, stdout=sps.PIPE)
    #     out = proc1.communicate()
    #     if proc1.returncode:
    #         logger.error("Failed: {}".format(out[0].decode()))
    #         raise OSError("Failed: {}".format(out[0].decode()))
    #     self.init_cache = out[0].decode().split(">")[-1].split("<")[0]

    #     # set new cache value to "true"
    #     if self.init_cache != "true":
    #         cmd0 = [
    #             self.binaries["vdb-config"], "--set", "cache-disabled:true"]
    #         out = sps.Popen(cmd0).communicate()
    #         logger.debug("vdb-config 'cache-disabled' set to 'true'")





    def binary_wrap(self):
        logger.warning("RUNNING")
        try:
            # self._set_vdbcfg()
            self._x_prefetch()
            self._x_fasterqd()
            self._x_kmerfreq()
            self._x_call_gce()
            # self.parse_results()

        finally:
            logger.info("removing tmp workdir")
            # self._set_vdbcfg(reset=True)
            self._clean_work()




if __name__ == "__main__":

    # SRS3758609    Ursus americanus    9643    11  SRR7811753
    SRR = "SRR7811753"
    tool = Kmunity(SRR, workdir="/home/deren/Downloads/kmunity-tmps")
    tool.binary_wrap()
