#!/home/GLBRCORG/mplace/anaconda3/bin/python
"""
@Program: rnaSeqCleanUp.py

@Purpose: Post Condor RNA-Seq pipeline run after all other condor jobs are
          complete in current directory.
          
@Input:   None 
          
                     
@Output: The following directories are created fastq, alignment, htseq, wig
         and the associated files are moved into them, i.e. bam files go 
         in the alignment directory.

@author: Mike Place
@Date:   3/1/2016
"""
import os
import sys

def cleanUp( cwd ):
    """
    Clean up and move output files.
    os.mkdir()
    os.rename( currentPath/, newPath )
    """
    cwd = cwd + "/"
    # move bam files to alignment directories
    os.mkdir( "alignment" )
    bamDir = cwd + "/alignment/"
    [ os.rename( (cwd + fn), (bamDir + fn) ) for fn in os.listdir(cwd) if fn.endswith(".bam") ]
    [ os.rename( (cwd + fn), (bamDir + fn) ) for fn in os.listdir(cwd) if fn.endswith(".bai") ] 
    # make wig directory
    os.mkdir( "wig" )
    widDir = cwd + "/wig/"
    [ os.rename( (cwd + fn), (widDir + fn) ) for fn in os.listdir(cwd) if fn.endswith(".wig.gz") ]
    # make HTseq directory
    os.mkdir( "htseq" )
    htsDir = cwd + "/htseq/"
    [ os.rename( (cwd + fn), (htsDir + fn) ) for fn in os.listdir(cwd) if fn.endswith("_HTseqOutput.txt") ]
    # make log file directory
    os.mkdir( "log" )
    logDir = cwd + "/log/"
    [ os.rename( (cwd + fn), (logDir + fn) ) for fn in os.listdir(cwd) if fn.endswith(".log") ]
    [ os.rename( (cwd + fn), (logDir + fn) ) for fn in os.listdir(cwd) if fn.endswith(".err") ]
    [ os.rename( (cwd + fn), (logDir + fn) ) for fn in os.listdir(cwd) if fn.endswith(".out") ]  
    # make a directory for sequence reads
    os.mkdir( "fastq" )
    seqDir = cwd + "/fastq/"
    [ os.rename( (cwd + fn), (seqDir + fn) ) for fn in os.listdir(cwd) if fn.endswith(".fastq") ]
    # delete sam files as they are no longer needed
    [ os.unlink(fn) for fn in os.listdir(cwd) if fn.endswith('.sam')]


def main():
    """
    Main 
    """
    currDir = os.getcwd()
    cleanUp(currDir)

if __name__ == "__main__":
    main()
 