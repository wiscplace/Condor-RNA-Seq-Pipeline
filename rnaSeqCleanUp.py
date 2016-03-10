#!/home/GLBRCORG/mplace/anaconda3/bin/python
"""
@Program: rnaSeqCleanUp.py

@Purpose: Post Condor RNA-Seq pipeline DAGman final process and clean up
          
@Input:  None           
                     
@Output: The following directories are created fastq, alignment, htseq, wig
         and the associated files are moved into them, i.e. bam files go 
         in the alignment directory. RPKM is calculated for all samples
         output is in RPKM.results.

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

def bam2wig():
    """
    Run bam2wig.pl, convert alignments from a Bam file into enumerated
    point data in a wig format. 
    bam2wig.pl --in run333.YPS163.10kreads.sort.bam --pos mid --strand --rpm --out YPS163.wig
    """
    program = '/opt/bifxapps/biotoolbox/scripts/bam2wig.pl'
    outFile = re.sub(r"bam", "wig", bamFile)
    cmd = [ program , '--in', bamFile, '--pos', 'mid', '--strand', '--rpm', '--out',  outFile]
    output = subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE ).communicate()
    result1 = output[0].decode( 'utf-8')
    result2 = output[1].decode( 'utf-8' )
    
    with open( 'bam2wig.log', 'a' ) as log:
        log.write("%s\n" %(outFile) )
        log.write(result1)
        log.write(result2)
        log.write("\n\n")    

def runRPKM( cwd, refer):
    """
    Run RPKM normalization on all HTSeq files in current working directory.
    
    Calls /home/GLBRCORG/mplace/scripts/RPKM.py

    RPKM.py -d <directory> -g REF.gff

    The assumption is that as part of this pipeline HTSeq has been run
    using CDS as the genome feature.
    
    """
    program = '/home/GLBRCORG/mplace/scripts/RPKM.py'
    gff     =  ref[refer][2] 
    cmd     = [ program, '-d', cwd, '-g', gff ]
    output  = subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE ).communicate()

def replaceWig( cwd ):
    """
    Run Kevin's script to replace chromosome names in wig file, for compatiblity with Mochi view.
    """
    wigDir  = cwd + "/wig/"
    os.chdir(wigDir)
    program = '/home/GLBRCORG/mplace/scripts/findreplace_WIG.pl'
    cmd     = [ program, wigDir ]
    output  = subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    result  = output[0].decode( 'utf-8' )
    with open( cwd + "/log/" + 'findreplace_WIG.log', 'w' ) as log:
        log.write(wigDir)
        log.write("\n")
        log.write("\n".join(cmd) )
        log.write("\n")
        log.write(result)


def main():
    """
    Main 
    """
    currDir = os.getcwd()
    runRPKM()
    bam2wig()
    replaceWig(currDir)
    cleanUp(currDir)

if __name__ == "__main__":
    main()
 
