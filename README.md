# reBLAST
Computational tool for identifying possible orthologues between genomes, given either DNA or Protein sequences.



#Requirements
- [Biopython](http://biopython.org/wiki/Main_Page)
- [NCBI BLAST](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)

#Installation
\- Checkout the source: `git clone https://github.com/kamanufred/reBLAST.git`  

#Usage
Run `python reBLAST.py -h` to see all the options.  

Usage: reBLAST.py [options]

    Options:
    -h, --help    show this help message and exit
    -a GENOME_1, --genome1=GENOME_1
                   First genome
    -b GENOME_2, --genome2=GENOME_2
                   Second genome
    -t MOL_TYPE, --type=MOL_TYPE
                   Molecule type [nucl or prot]
    -o OUTPUT_FILE, --output=OUTPUT_FILE
                  Output file with orthologs

#Contact
Contact `frederick(dot)kamanu(at)gmail.com` for feedback regarding the software.  

#License    
This software is provided under the GNU General Public License. See the included file.