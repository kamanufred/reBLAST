#!/usr/bin/env python

# Copyright 2015 Frederick Kinyua Kamanu, Ph.D [frederick dot kamanu at gmail.com]
# All Rights Reserved.

#*****************************************************************************************
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>.
#*****************************************************************************************

import sys
import os
import subprocess
import shutil
from optparse import OptionParser
from Bio import SeqIO
from Bio.Blast import NCBIXML


#ORTHOgrab main class
#-----------------------------------------------------------------------------------------
class OrthoGrab(object):
    """
    Run a blast analysis between any two genomes as input and output
    a list of orthologues genes
    """
    def __init__(self,genome1,genome2,moltype,results):
        """
        Initialize the class
        """
        self.genome1 = genome1
        self.genome2 = genome2
        self.moltype = moltype
        self.results = results


    def CreateBlastDir(self):
        """
        Create a directory for blastdb directory and copy data file
        """
        temp_dir = 'blastdb'
        try:
            os.mkdir(temp_dir)
        except OSError:
            shutil.rmtree(temp_dir)
            os.mkdir(temp_dir)


    def CreateBlastDb(self):
        """
        Build a blast database
        """
        blastdb_a = 'blastdb/a'
        blastdb_b = 'blastdb/b'
        logfile = 'temp.log'
        seq = ''
        if self.moltype == 'prot':
            seq = 'prot'
        elif self.moltype == 'nucl':
            seq = 'nucl'
        else:
            sys.exit("\nPlease enter 'nucl' or 'prot' for molecule type\n")

        os.system("makeblastdb -in %s -dbtype %s -out %s >logfile" % (self.genome1,seq,blastdb_a))
        os.system("makeblastdb -in %s -dbtype %s -out %s >logfile" % (self.genome2,seq,blastdb_b))


    def RunBlast(self):
        """"
        execute blast
        """
        program = ''
        query1 = self.genome1
        query2 = self.genome2
        database1 = 'blastdb/a'
        database2 = 'blastdb/b'
        outfile1 = 'blastout1.xml'
        outfile2 = 'blastout2.xml'
        outformat = 5
        evalue = 0.001
        max_target_seqs = 5
        num_threads = 2

        if self.moltype == 'nucl':
            program = 'blastn'
        else:
            program = 'blastp'

        sys.stderr.write("\nRunning a vs b %s ...\n" % (program))
        os.system("%s -query %s -db %s -out %s -outfmt %d -evalue %f \
        -max_target_seqs %d -num_threads %d" % (program,query1,database2,\
        outfile1,outformat, evalue,max_target_seqs,num_threads))

        sys.stderr.write("\nRunning b vs a %s ...\n" % (program))
        os.system("%s -query %s -db %s -out %s -outfmt %d -evalue %f \
        -max_target_seqs %d -num_threads %d" % (program,query2,database1,\
        outfile2,outformat, evalue,max_target_seqs,num_threads))


    def __BlastXMLParser(self,blast_file):
        """
        Parse BlastXML results and return a dictionary of query <=> tophit
        """
        matches = {}
        result_handle = open(blast_file,'rU')
        blast_records = NCBIXML.parse(result_handle)
        for blast_record in blast_records:
            query = blast_record.query
            if blast_record.alignments:
                tophit_id = blast_record.alignments[0].title.split(" ")[0]
                tophit_def = " ".join(blast_record.alignments[0].title.split(" ")[1:])
                query_id = query.split("|")[1]
                tophit_id_clean = tophit_def.split("|")[1]
                matches[query_id] = tophit_id_clean
        return matches


    def GetOrthologs(self):
        """
        Parse the blast results
        """
        sys.stderr.write("\nParsing BlastXML results ...\n")
        matches1 = self.__BlastXMLParser('blastout1.xml')
        matches2 = self.__BlastXMLParser('blastout2.xml')

        write_handle = open(self.results,"w")
        for key,value in matches1.items():
            try:
                if matches2[value] == key:
                    write_handle.write("%s\t%s\n" % (key,value))
            except KeyError:
                pass


    def CleanUp(self):
        """
        Simple clean up all the temp files
        """
        os.system("rm -rf blastdb blastout1.xml blastout2.xml logfile")
        sys.stderr.write("\nFinal results written to %s\n\n" % (self.results))



#Test for dependency installation
#-----------------------------------------------------------------------------------------
def Which(program):
    """
    Check if a program has been installed in a *NIX system and is in the path
    """
    status = 0
    try:
    # pipe output to /dev/null for silence
        null = open("/dev/null", "w")
        subprocess.Popen(program, stdout=null, stderr=null)
        null.close()
        status = 1

    except OSError:
        status = 0

    return status


#Parsing commandline options
#-----------------------------------------------------------------------------------------
def CommandlineOptions():

    parser = OptionParser(usage="usage: %prog [options]")
    parser.add_option("-a", "--genome1",
                      type="string",
                      dest="genome_1",
                      help="First genome")
    parser.add_option("-b", "--genome2",
                      type="string",
                      dest="genome_2",
                      help="Second genome")
    parser.add_option("-t", "--type",
                      type="string",
                      dest="mol_type",
                      help="Molecule type [nucl or prot]")
    parser.add_option("-o", "--output",
                      type="string",
                      dest="output_file",
                      help="Output file with orthologs")

    (options, args) = parser.parse_args()
    options_args_parser = [options,args,parser]

    return options_args_parser

#Main
#-----------------------------------------------------------------------------------------
def main():

    options, args, parser = CommandlineOptions()
    genome1 = options.genome_1
    genome2 = options.genome_2
    moltype = options.mol_type
    output = options.output_file

    if genome1 is None or genome2 is None or moltype is None or output is None:
        print "\nError: A mandatory option is missing!\n"
        parser.print_help()
        sys.exit(-1)

    else:
        if Which('blastn') == 0:
            sys.exit("\nBlast not installed!\n")
        else:
            OrthoGrabObject = OrthoGrab(genome1,genome2,moltype,output)
            OrthoGrabObject.CreateBlastDir()
            OrthoGrabObject.CreateBlastDb()
            OrthoGrabObject.RunBlast()
            OrthoGrabObject.GetOrthologs()
            OrthoGrabObject.CleanUp()


#Execute main
#-----------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()