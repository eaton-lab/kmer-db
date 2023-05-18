#!/usr/bin/env python

""" command line tool for kmunity"""

from pkg_resources import get_distribution
import argparse
import sys
import kmunity


EPILOG = """
  * Select a Run ID and write results to mammals dir.
    kmunity -s SRR7811753 -d mammals -w /scratch/ -r ./kmunity

  * Autoselect a Run ID to contribute to mammals dir.
    kmunity -d mammals -w /scratch/ -r ./kmunity

  * Documentation: http://github.com/eaton-lab/kmunity
    """


class CLI:
    """Command line tool for Kmunity GCE analysis.
    """
    def __init__(self):

        # the parser object will be used to fill args, parsedict, and load data
        self.args = None
        self.parsedict = None
        self.data = None

        # parse the command line arguments
        self.parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog=EPILOG)
        self._parse_command_line()
        self.args = self.parser.parse_args()

        # configure kmunity and sratools
        if self.args.config:
            self.config()

        # finally run the requested functions
        else:
            self.run()

    def _parse_command_line(self):
        """ Parse CLI args."""

        # if no args then return help message
        if len(sys.argv) == 1:
            self.parser.print_help()

        # add arguments
        self.parser.add_argument(
            '-v', '--version',
            action='version',
            version=str(get_distribution('kmunity')),
        )

        self.parser.add_argument(
            "-s", dest="srr", type=str, default=None,
            help="SRR Run ID")

        self.parser.add_argument(
            '-r', dest="repo", type=str, default=".",
            help="path to forked kmunity repository.")

        self.parser.add_argument(
            '-w', dest="workdir", type=str, default=".",
            help="path to working directory where tmp files will be stored.")

        self.parser.add_argument(
            "-d", dest="database", type=str,
            default="mammals",
            help="database to contribute to (mammals, birds, ...)")

        self.parser.add_argument(
            "--config", action='store_true',
            help="configure kmunity with sra-tools")

        self.parser.add_argument(
            "--search", action='store_true',
            help="search for an SRR from db")

        # self.parser.add_argument(
        #     "-t", dest="tmpdir", type=str,
        #     default="/tmp/",
        #     help="path to a director where tmp software will be stored."
        #     )

    def config(self):
        kmunity.Kmunity(config=True)

    def run(self):
        tool = kmunity.Kmunity(
            self.args.srr,
            self.args.database,
            self.args.workdir,
            self.args.repo,
        )
        if not self.args.search:
            tool.binary_wrap()


def main():
    CLI()


if __name__ == "__main__":
    main()
