#!/home/GLBRCORG/mplace/anaconda3/bin/python
"""
@Program: RPKM.py

@Purpose: Run RPKM normalization on HTSeq results.

@Input:   HTSeq file, Reference GFF, gff feature type i.e. "CDS", "gene", genome name
          default feature type is CDS, as the RNA-Seq pipeline runs HTSeq with CDS.
          
          HTSeq file -- just standard output from HTSeq
              ( see: http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html )
              
          Reference GFF  -- General Feature Format, usually comes with Reference Genome
          
          gff feature type:  default is CDS, but gene and others may be used, 
                             see your gff file for what is available.
          genome name: R64-1-1 -- SGD S. cerevisiae genome R64-1-1
                       R64-2-1 -- SGD S. cerevisiae genome R64-2-1
                       PAN     -- PanGenome
                       Y22-3   -- GLBRC Y22-3 S. cerevisiae genome
                
@Output:  RPKM results for a feature type in a tab delimited table.
          First column is the Gene (feature) Name
          the next columns are a strains RPKM value for a feature.
          
          If a gene is not found in the GFF file, the RPKM is set to zero.

example:

Gene    Gasch136_test.sort.bam_HTseqOutput.txt  Gasch138_test.sort.bam_HTseqOutput.txt
Q0010   0.0     0.0
Q0017   0.0     0.0
Q0032   0.0     0.0
Q0045   2.009926923942714       6.312318612212933
Q0050   0.397603231814493       0.45727088358354523
Q0055   0.0     0.22849861627304624
Q0060   0.0     0.0
         
@Dependencies: 
        Python 3, 
        GFFUtils 
        Samtools 
      
@Output: Tab delimited text table, columns are strains, rows are gene names
@author: Mike Place
@Date:   10/1/2015

"""
import argparse	
import math
import os
import re
import subprocess
import sys
from collections import defaultdict

featureDict     = {}            # key = gene name,   value = tuple(Gene_LENGTH, StartPos, StopPos)
totalReadCounts = {}            # key = sample name, value = total aligned read counts
bamList         = []            # list of sample bam files to process
genomeSize      = 0
geneList        = []            # list of genes in HTSeq files, order = HTSeq file order
result          = defaultdict(dict)  # dict of dicts key = strain , second key = gene name value =  RPKM value

def getGeneList(file):
    """
    Get the list of genes in the same order as the HTSeq files.
    for outputing results in HTSeq order
    """
    htseqName = file + "_HTseqOutput.txt"
    
    with open(htseqName) as htseq:
        for line in htseq:
            if line.startswith('_'):
                continue
            row = line.split()
            geneList.append(row[0])

def getTotalReadCounds( bamList  ):
    """
    bamList = list of all files names .bam in current working directory.
    Open and sum the read counts in each htseq file.       
    This will include all read counts for each gene as well as the counts
    for the non-gene categories at the bottom of column:
    __no_feature    
    __ambiguous     
    __too_low_aQual 
    __not_aligned   
    __alignment_not_unique  
    
    """
    for b in bamList:
        htseqName = b + "_HTseqOutput.txt"
        
        # calculate the total number of aligned reads        
        totalAlignedReads = 0
        
        if not checkFormat(htseqName):
            with open('RPKM.log', 'a') as log:
                log.write("\n\nHTSeq file format ERROR:\t")
                log.write("%s \n" % htseqName)
                log.write("\tHTSeq format expects 2 columns: gene_name, gene_count\n")
                log.write("\tMoving on to next file.")
                continue        
        # sum read counts
        with open(htseqName,'r') as htseq:
            for x in htseq:
                x = x.strip()               # get rid of that pesky newline
                row = x.split('\t')               
                totalAlignedReads += float(row[1])
        
        totalReadCounts[htseqName] = totalAlignedReads
        
    with open('RPKM.log','a') as log:
        for strain,count in totalReadCounts.items():
            log.write("\n")
            log.write("%s has %d reads aligned to genes." %(strain, count))
    log.close()
    
def createGeneDict( gffFile, feature ):
    """
    Create a dictionary of genes and gene lengths 
    key = gene_name value = tuple of (length)
    Used in the RPKM Calculation
    """            
    #Open GFF file
    with open(gffFile, 'r') as gff:
        for g in gff:
            if g.startswith('#'):            # skip comment
                continue
            if g.startswith('>'):            # if fasta at end of gff lines, stop loop, assumes fasta after annotations
                break
            grow = g.split()
            if grow[2] == feature:            # only process lines with gene feature type
                info = grow[8].split(';')
                for i in info:
                    if i.startswith('Name'):
                        geneName = i.split('=')
                        geneLen  = float(int(grow[4]) - int(grow[3])) + 1.0
                        # a feature is maybe listed more than once, so sum all of them
                        if geneName[1] in featureDict:
                            featureDict[geneName[1]] += (geneLen)
                        else:
                            featureDict[geneName[1]] = (geneLen)                                
                        
    with open('RPKM.log','a') as log:
        dictLength = len(featureDict)
        log.write("\nNumber of Genes in GFF: %d\n" % dictLength)
        log.close()
    
def calcRPKM(genomeSize):
    """
    use HTSeq results files for each bam to calculate RPKM.
    
    RPKM(x) = (10^9 * C)/(N * L)
    
    C = Number of reads mapped to gene x
    N = Total number of aligned reads
    L = Length of gene
    """    
    
    for file, count in totalReadCounts.items():
        if not os.path.exists(file):
            print("\n\tHTSeq file does not exist.\n")
            sys.exit(1)
        with open(file,'r') as htseq:
            for line in htseq:
                if line.startswith('_'):
                    continue
                line = line.rstrip()
                gene = line.split()         # gene[0] == geneName  gene[1] == HTSeq alignment count
                
                if gene[0] in featureDict:  # featureDict[geneName] == (GeneLength)
                    if count == 0 or featureDict[gene[0]] == 0 or gene[1] == 0:
                        rpkm = 0
                    else:
                        rpkm = (float(math.pow(10,9)) * float(gene[1]))/(float(featureDict[gene[0]]) * float(count))
                    if not result[file]:
                        result[file][gene[0]] = rpkm
                    else:
                        result[file][gene[0]] = rpkm
                else:
                    result[file][gene[0]] = 0

def getRefGenomeSize( gen ):
    """
    Get reference genome size, values computed from the fasta files
    for each genome.    
    """
    
    R64-1-1 = 12157105
    R64-2-1 = 12157105
    Pan     = 13345927
    Y22     = 12170480

    
    genomeSize = 0
    
    
    return genomeSize
                    
def checkFormat(htseqFile):
    """
    Make sure the HTSeq file is in the expected format.
    Just a dumb check looking for 2 columns.
    """
    with open(htseqFile, 'r') as file:
        line = file.readline()
        test = line.split()
        if len(test) == 2:
            return True
        else:
            return False

def main():
    """
    Main 
    """
    # Variables
    htseqFiles   = []         # cds files to process, in curr working dir
    
    cmdparser = argparse.ArgumentParser(description="Run RPKM on all HTSeq files in Current Directory.",
                                        usage='%(prog)s -d <directory> -g REF.gff [-f gene] ' ,prog='RPKM.py'  )                                  
    cmdparser.add_argument('-d', '--dir',  action='store',      dest='DIR' , help='Directory path containing the HTSeq output files.', metavar='')
    cmdparser.add_argument('-g', '--gff',  action='store',      dest='GFF' , help='GFF file to use.', metavar=''  )
    cmdparser.add_argument('-f', '--feature', action='store',   dest='FEATURE', help='Feature type to use w/ GFF file, default = CDS.', metavar='')        
    cmdparser.add_argument('-i', '--info', action='store_true', dest='INFO', help='Detailed description of program.')
    cmdResults = vars(cmdparser.parse_args())
    
    # if no args print help
    if len(sys.argv) == 1:
        print("")
        cmdparser.print_help()
        sys.exit(1)
             
    if cmdResults['INFO']:
        print("\n  RPKM.py -d /Full/Path/To/HTSeq_output -g Reference GFF file [-f feature]")
        print("\n  Purpose: Run RPKM normalization on HTSeq read counts.")
        print("\n  Input  : Full path to HTSeq results directory.")
        print("\n  The bam and bam index files need to be in current directory for script to work.")
        print("\n  Output : Single RPKM results file, where each column denotes a strain.")    
        print("\n  Usage  : RPKM.py -d /home/GLBRCORG/user/My_HTSeq_Dir -g REF.gff [-f feature] ")        
        print("  ")       
        print("\tTo see Python Docs and get a better explaination of the program:")
        print("\n\tOpen python console and enter")
        print("\timport sys")
        print("\tsys.path.append('/full/path/to/script')")
        print("\timport RPKM")
        print("\thelp(RPKM)")
        print("\n\tSee Mike Place for any problems or suggestions.\n")
        sys.exit(1)
        
    if cmdResults['DIR']:
        htSeqDir = cmdResults['DIR']
    
    if not os.path.exists(htSeqDir):
        print("\n\t-d HTSeq directory does not exist.\n")
        cmdparser.print_help()
        sys.exit(1)
    else:
        # Here we are assuming the .fai index file for reference exists
        for f in os.listdir(htSeqDir):
            if re.search(r'\_HTseqOutput.txt$',f):
                htseqFiles.append(f)
    
    if cmdResults['GFF']:
        gffFile = cmdResults['GFF']
        if not os.path.exists(gffFile):
            print("\n\t -g GFF file does not exist.\n")
            cmdparser.print_help()
            sys.exit(1)
    else:
        print("\n\t -g GFF file does not exist.\n")
        cmdparser.print_help()
        sys.exit(1)
        
    
    if cmdResults['FEATURE']:
        feature = cmdResults['FEATURE']
    else:
        feature = 'CDS'    
    
    # Log the initial parameters
    with open('RPKM.log','w') as log:
        log.write("RPKM running in the following Directory: %s\n" %(htSeqDir))
        log.write("Processing the following files:\n")
        for ht in htseqFiles:
            log.write("\t\t%s\n" %(ht))
        log.write("\nGFF file     : %s\n" %(gffFile))
        log.write("\nUsing Feature: %s\n" %(feature))
        log.close()
        
    # Get a list of bam files in current directory
    for bam in os.listdir():
        if bam.endswith('.bam'):
            bamList.append(bam)
        
    # if there are no files with the .bam extension, print help and exit.
    if not bamList:
        print("No bam files found in current directory.")
        print("Exiting program")
        cmdparser.print_help()
        sys.exit()  
    
    # Start Real Work
    createGeneDict(gffFile, feature)
    getTotalReadCounds(bamList)
    genomeSize = getRefGenomeSize(bamList[0])  
    getGeneList(bamList[0])
    calcRPKM(genomeSize)
    
    with open('RPKM.results','w') as out:
    #print header
        header = bamList
        header.insert(0,"Gene")
        newHeader = "\t".join(header)
        out.write(newHeader + "\n")        
        
        for g in geneList:
            data = []
            data.append(g)
            out.write("%s\t" %(g))
            row = []
            for b in bamList:
                name = b + "_HTseqOutput.txt"
                if name in result:
                    row.append(str(result[name][g]))
            newRow = "\t".join(row)
            out.write(newRow + "\n")
                  


if __name__ == "__main__":
    main()
    




