#!/home/GLBRCORG/mplace/anaconda3/bin/python
"""
@Program: RPKM.py

@Purpose: Run RPKM normalization on HTSeq results on the HTCondor RNA-Seq pipeline.
          Simplified command line.

@Input:   HTSeq file, gff feature type i.e. "CDS", "gene", genome name
          default feature type is CDS, as the RNA-Seq pipeline runs HTSeq with CDS.
          
          HTSeq file -- just standard output from HTSeq
              ( see: http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html )
                      
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
htList         = []            # list of sample bam files to process
genomeSize      = 0
geneList        = []            # list of genes in HTSeq files, order = HTSeq file order
result          = defaultdict(dict)  # dict of dicts key = strain , second key = gene name value =  RPKM value

# list of available gff files
refGFF = { 'R64' : "/home/GLBRCORG/mplace/data/reference/S288C_reference_genome_R64-1-1_20110203/saccharomyces_cerevisiae_R64-1-1_20110208_noFasta.gff",
           'Y22' : "/home/GLBRCORG/mplace/data/reference/Y22-3/Y22-3_Final_GFF.gff",
           'PAN' : "/home/GLBRCORG/mplace/data/reference/PanGenome/PanGenome-Final-R64-2-1.gff",
           'R64-2' : "/home/GLBRCORG/mplace/data/reference/S288C_reference_genome_R64-2-1_20150113/saccharomyces_cerevisiae_R64-2-1_20150113_noFasta.gff"                    
    }

def getGeneList(file):
    """
    Get the list of genes in the same order as the HTSeq files.
    for outputing results in HTSeq order
    """
    htseqName = file
    
    with open(htseqName) as htseq:
        for line in htseq:
            if line.startswith('_'):
                continue
            row = line.split()
            geneList.append(row[0])

def getTotalReadCounds( htList  ):
    """
    htList = list of all HTSeq files in current working directory.
    Open and sum the read counts in each htseq file.       
    This will include all read counts for each gene as well as the counts
    for the non-gene categories at the bottom of column:
    __no_feature    
    __ambiguous     
    __too_low_aQual 
    __not_aligned   
    __alignment_not_unique  
    
    """
    for ht in htList:
        htseqName = ht
        
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
    genomeDict = { 'R64' : 12157105, 'R64-2' : 12157105, 'PAN' : 13345927, 'Y22' : 12170480 }
    return genomeDict[gen]
                    
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
    cmdparser.add_argument('-f', '--feature', action='store',   dest='FEATURE', help='Feature type to use w/ GFF file, default = CDS.', metavar='')
    cmdparser.add_argument('-r', '--reference', action='store', dest='REFERENCE', help='Genome Reference for genome size information.', metavar='')
    cmdparser.add_argument('-i', '--info', action='store_true', dest='INFO', help='Detailed description of program.')
    cmdResults = vars(cmdparser.parse_args())
    
    # if no args print help
    if len(sys.argv) == 1:
        print("")
        cmdparser.print_help()
        sys.exit(1)
             
    if cmdResults['INFO']:
        print("\n  RPKM.py -r genome [-f feature] ")
        print("\n  Purpose: Run RPKM normalization on HTSeq read counts.")
        print("\n  Input  : reference genome name.")
        print("\n  Output : Single RPKM results file, where each column denotes a strain.")    
        print("\n  Usage  : RPKM.py -r genome [-f feature] ")
        print("\n           genome options R64 (SGD R64-1-1), R64-2 (SGD R64-2-1), Y22-3 (GLBRC), PAN (PanGenome)")
        print("  ")       
        print("\tTo see Python Docs and get a better explaination of the program:")
        print("\n\tOpen python console and enter")
        print("\timport sys")
        print("\tsys.path.append('/full/path/to/script')")
        print("\timport RPKM")
        print("\thelp(RPKM)")
        print("\n\tSee Mike Place for any problems or suggestions.\n")
        sys.exit(1)
        
    htSeqDir = os.getcwd()
    # Here we are assuming the .fai index file for reference exists
    for f in os.listdir(htSeqDir):
        if re.search(r'\_HTseqOutput.txt$',f):
            htseqFiles.append(f)
        
    if cmdResults['FEATURE']:
        feature = cmdResults['FEATURE']
    else:
        feature = 'CDS'

    if cmdResults['REFERENCE'] == 'R64-2' or cmdResults['REFERENCE'] == 'PAN' or cmdResults['REFERENCE'] == 'Y22':
        refGenome = cmdResults['REFERENCE']
        gffFile  = refGFF[refGenome]
    else:
        refGenome = 'R64'
        gffFile = refGFF['R64']
    
    # Log the initial parameters
    with open('RPKM.log','w') as log:
        log.write("RPKM running in the following Directory: %s\n" %(htSeqDir))
        log.write("Processing the following files:\n")
        for ht in htseqFiles:
            log.write("\t\t%s\n" %(ht))
        log.write("\nGFF file     : %s\n" %(gffFile))
        log.write("\nUsing Feature: %s\n" %(feature))
        log.write("\nUsing Genome : %s\n" %(refGenome))
        log.close()

    for ht in os.listdir():
        if ht.endswith('_HTseqOutput.txt'):
            htList.append(ht)

    if not htList:
        print("No HTSeq files found in current directory")
        print("Exiting Program")
        
    # Start Real Work  
    createGeneDict(gffFile, feature)
    getTotalReadCounds(htList)
    genomeSize = getRefGenomeSize(refGenome)
    getGeneList(htList[0])
    calcRPKM(genomeSize)
    
    with open('RPKM.results','w') as out:
    #print header
        header = htList
        header.insert(0,"Gene")
        newHeader = "\t".join(header)
        out.write(newHeader + "\n")        
        
        for g in geneList:
            data = []
            data.append(g)
            out.write("%s\t" %(g))
            row = []
            for ht in htList:
                name = ht 
                if name in result:
                    row.append(str(result[name][g]))
            newRow = "\t".join(row)
            out.write(newRow + "\n")

if __name__ == "__main__":
    main()
    




