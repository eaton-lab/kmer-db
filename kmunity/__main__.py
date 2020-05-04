#!/usr/bin/env python

""" command line tool for kmunity"""


from __future__ import print_function

from pkg_resources import get_distribution
import argparse
import sys




class CLI:
    def __init__(self):


        # the parser object will be used to fill args, parsedict, and load data
        self.args = None
        self.parsedict = None
        self.data = None
        self.parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog=EPILOG)


        # parse the command line to get CLI args
        self._parse_command_line()
        self.args = self.parser.parse_args()

        # bail if no args for 'params' or 'new'
        self.check_args()

        # check for download argument
        if self.args.download:
            self._flagdownload()
            sys.exit(0)

        # check for merge of branches
        if self.args.merge:
            self.merge_assemblies()
            sys.exit(0)

        # finally run the requested functions
        self.run()


    def _parse_command_line(self):
        """ Parse CLI args."""

        # if no args then return help message
        if len(sys.argv) == 1:
            self.parser.print_help()

        ## add arguments 
        self.parser.add_argument(
            '-v', '--version', 
            action='version', 
            version=str(get_distribution('ipyrad')),
        )
        self.parser.add_argument('-r', "--results", action='store_true',
            help="show results summary for Assembly in params.txt and exit")

        self.parser.add_argument('-f', "--force", action='store_true',
            help="force overwrite of existing data")

        self.parser.add_argument('-q', "--quiet", action='store_true',
            help="do not print to stderror or stdout.")

        self.parser.add_argument('-d', "--debug", action='store_true',
            help="print lots more info to ipyrad_log.txt.")

        self.parser.add_argument('-n', dest="new", type=str, default=None, 
            help="create new file 'params-{new}.txt' in current directory")

        self.parser.add_argument('-p', dest="params", type=str, default=None,
            help="path to params file for Assembly: params-{assembly_name}.txt")

        self.parser.add_argument('-s', dest="steps", type=str, default=None,
            help="Set of assembly steps to run, e.g., -s 123")

        self.parser.add_argument('-b', dest="branch", type=str, default=None, 
            nargs="*",
            help="create new branch of Assembly as params-{branch}.txt, and " + \
            "can be used to drop samples from Assembly.")

        self.parser.add_argument('-m', dest="merge", default=None, nargs="*",
            help="merge multiple Assemblies into one joint Assembly, and " + \
            "can be used to merge Samples into one Sample.")

        self.parser.add_argument("-c", metavar="cores", dest="cores",
            type=int, default=0,
            help="number of CPU cores to use (Default=0=All)")

        self.parser.add_argument("-t", metavar="threading", dest="threads",
            type=int, default=2,
            help="tune threading of multi-threaded binaries (Default=2)")


        self.parser.add_argument("--download", dest="download", type=str, 
        nargs="*", default=None,  # const="default",
            help="download fastq files by accession (e.g., SRP or SRR)")




HEADER = """
 -------------------------------------------------------------
  ipyrad [v.{}]
  Interactive assembly and analysis of RAD-seq data
 -------------------------------------------------------------\
 """.format(str(get_distribution('ipyrad')).split()[1])
# ip.__version__)


EPILOG = """
  * Example command-line usage: 
    ipyrad -n data                       ## create new file called params-data.txt 
    ipyrad -p params-data.txt -s 123     ## run only steps 1-3 of assembly.
    ipyrad -p params-data.txt -s 3 -f    ## run step 3, overwrite existing data.

  * HPC parallelization across 32 cores
    ipyrad -p params-data.txt -s 3 -c 32 --MPI

  * Print results summary 
    ipyrad -p params-data.txt -r 

  * Branch/Merging Assemblies
    ipyrad -p params-data.txt -b newdata  
    ipyrad -m newdata params-1.txt params-2.txt [params-3.txt, ...]

  * Subsample taxa during branching
    ipyrad -p params-data.txt -b newdata taxaKeepList.txt

  * Download sequence data from SRA into directory 'sra-fastqs/' 
    ipyrad --download SRP021469 sra-fastqs/ 

  * Documentation: http://ipyrad.readthedocs.io
    """



def main():
    CLI()


if __name__ == "__main__": 
    main()
