#!/usr/bin/env python

"""
Download fastq files and calculate heterozygosity with gce
"""


import os
import sys
import subprocess as sps



class KmerHet:
    """
    By default it looks for software in your PATH. 
    """
    def __init__(self, srr, workdir, **kwargs):

        # store args
        self.srr = srr

        # check kwargs: e.g., user-supplied binary paths
        self.kwargs = {
            "prefetch": None,
            "fasterq-dump": None,
            "kmerfreq": None,
            "gce": None,
            "outdir": os.path.realpath(os.path.expanduser(workdir))
        }
        self._check_kwargs()

        # check binaries/paths
        self._check_binaries()


    def _check_binaries(self):
        for binary in ["prefetch", "fasterq-dump", "kmerfreq", "gce"]:
            if not self.kwargs[binary]:
                self.kwargs[binary] = os.path.join(sys.prefix, "bin", binary)
            if not os.path.exists(self.kwargs[binary]):
                raise ImportError("missing dependency: {}".format(binary))


    def _check_kwargs(self):
        if not os.path.exists(self.kwargs['outdir']):
            os.makedirs(self.kwargs['outdir'])


    def call_prefetch(self):
        # print the version
        proc = sps.Popen(
            [self.kwargs["prefetch"], "-V"], 
            stderr=sps.STDOUT, 
            stdout=sps.PIPE,
        )
        out = proc.communicate()[0].split()[-1]
        print("[srrhet] prefetch {}".format(out))

        # run prefetch
        cmd = [
            self.kwargs["prefetch"], self.srr, 
            "-O", self.kwargs["outdir"],
        ]
        proc = sps.Popen(cmd, stderr=sps.STDOUT, stdout=sps.PIPE)
        out = proc.communicate()
        if proc.returncode:
            raise Exception(out[0])


    def call_fasterq_dump(self):
        # print the version
        proc = sps.Popen(
            [self.kwargs["fasterq-dump"], "-V"],
            stderr=sps.STDOUT, 
            stdout=sps.PIPE,
        )
        out = proc.communicate()[0].split()[-1]
        print("[srrhet] fasterq-dump {}".format(out))

        # call the tool
        cmd = [
            self.kwargs["fasterq-dump"], self.srr, 
            "-O", self.kwargs["outdir"],
        ]
        proc = sps.Popen(cmd, stderr=sps.STDOUT, stdout=sps.PIPE)
        out = proc.communicate()
        if proc.returncode:
            raise Exception(out[0])

        # write a tmp SRR.lib file
        libfile = os.path.join(self.workdir, "{}_files.lib".format(self.srr))
        with open(libfile, 'w') as out:
            out.write("...")


    def call_kmerfreq(self):
        proc = sps.Popen(
            [self.kwargs["kmerfreq"], "-h"],
            stderr=sps.STDOUT, 
            stdout=sps.PIPE,
        )
        out = proc.communicate()[0].decode().split("\n")
        out = [i for i in out if "Version" in i]
        vers = out.split()[-1]
        print("[srrhet] fasterq-dump {}".format(vers))

        # call the tool
        cmd = [
            self.kwargs["kmerfreq"],
            "-k", "17",
            "-t", "4",
            "-p", self.srr,
            "{}_files.lib".format(self.srr),
        ]
        proc = sps.Popen(cmd, stderr=sps.STDOUT, stdout=sps.PIPE)
        out = proc.communicate()
        if proc.returncode:
            raise Exception(out[0])


    def call_gce(self):
        # print the version
        proc = sps.Popen(
            [self.kwargs["fasterq-dump"], "-V"],
            stderr=sps.STDOUT, 
            stdout=sps.PIPE,
        )
        out = proc.communicate()[0].split()[-1]
        print("[srrhet] fasterq-dump {}".format(out))


    def parse_results(self):
        pass


    def run(self):
        self.call_prefetch()
        self.call_fasterq_dump()
        self.call_kmerfreq()
        self.call_gce()
        self.parse_results()




if __name__ == "__main__":

    # SRS3758609    Ursus americanus    9643    11  SRR7811753
    SRR = "SRR7811753"
    tool = HetSRR(SRR, outdir="/home/deren/Downloads/SRRHET")
    tool.run()
