#!/home/GLBRCORG/mplace/anaconda3/bin/python
"""
@Program: HTseqReport.py

@Purpose: Create a simple report for each HTseq result file.

@Input:   HTSeq files 
                          
@Output: Tab delimited text table
filename    Total_Reads    Num_Reads_Aligned_to_Feature    %Aligned_to_Feature
    
@author: Mike Place
@Date:   5/1/2016
"""
import argparse	
import os
import re
import sys
from collections import defaultdict

counts = defaultdict(list)

def getReadCounts( counts, htList  ):
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

    Also sum only the aligned read counts
    return Total reads and Aligned read count
    """
    for ht in htList:
        htseqName = ht
        
        # calculate the total number of aligned reads        
        totalReads   = 0
        alignedReads = 0
              
        # sum read counts
        with open(htseqName,'r') as htseq:
            for x in htseq:
                x = x.strip()               # get rid of that pesky newline
                row = x.split('\t')
                totalReads += float(row[1])
                if x.startswith('__'):
                    continue
                else:
                    alignedReads += float(row[1])
                    
        percentAligned = (alignedReads/totalReads) * 100
        counts[ht].append(totalReads)
        counts[ht].append(alignedReads)
        counts[ht].append(percentAligned)

def main():
    """
    Main 
    """
    # Variables
    htseqFiles   = []         # htseq files to process
    
    cmdparser = argparse.ArgumentParser(description="Produce report on all HTSeq files in current directory.",
                                        usage='%(prog)s ' ,prog='HTseqReport.py'  )                                  
    cmdparser.add_argument('-i', '--info', action='store_true', dest='INFO', help='Detailed description of program.')
    cmdResults = vars(cmdparser.parse_args())
             
    if cmdResults['INFO']:
        print("\n  HTseqReport.py ")
        print("\n  Purpose: Create a simple report of HTSeq results.")
        print("\n  Input  : none, just need to be in directory w/ HTseq results files")
        print("           Files expected to end with _HTseqOutput.txt.gz ")
        print("\n  Output : Single tab delimited report file, 4 colums, filename    Total_Reads    Num_Reads_Aligned_to_Feature    %Aligned_to_Feature")    
        print("\n  Usage  : HTseqReport.py ")
        print("  ")       
        print("\tTo see Python Docs and get a better explaination of the program:")
        print("\n\tOpen python console and enter")
        print("\timport sys")
        print("\tsys.path.append('/full/path/to/script')")
        print("\timport HTseqReport")
        print("\thelp(HTseqReport)")
        print("\n\tSee Mike Place for any problems or suggestions.\n")
        sys.exit(1)
        
    htSeqDir = os.getcwd()
    for f in os.listdir(htSeqDir):
        if re.search(r'\_HTseqOutput.txt.gz$',f):
            htseqFiles.append(f)

    if not htseqFiles:
        print("No HTSeq files found in current directory")
        print("Exiting Program")
        sys.exit(1)

    getReadCounts(counts, htseqFiles)
    print(counts)
    
    with open('HTSeq-report.txt','w') as out:
    #print header
        header = ['filename', 'Total_Reads', 'Reads_Aligned_to_Feature', '%Aligned_to_Feature' ]
        Header = "\t".join(header)
        out.write(Header + "\n")        
        
        for k,v in counts.items():
            row = "\t".join([str(x) for x in v])
            out.write("%s\t" %(k))
            out.write("%s\n" %(row))   

if __name__ == "__main__":
    main()
    




