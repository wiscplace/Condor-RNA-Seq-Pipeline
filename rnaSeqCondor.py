#!/home/GLBRCORG/mplace/anaconda3/bin/python
"""
Program: rnaSeqCondor.py

Purpose: Implement the currently used Gasch lab RNA-Seq pipeline using condor
         dagman for job managment. This is a reimplementation of the original
         RNA-Seq pipeline (https://gitlab.wei.wisc.edu/mplace/RNA-Seq-Pipeline)

         HTCondor https://research.cs.wisc.edu/htcondor/

Input:  reference file, single or paired

    text file with a list of RNA-Seq fastq files to be processed
    one file name per line. 
    to generate:  /bin/ls *.fastq > input.txt
    optional parameters: -htseq reverse   ( for HTSeq )
   
Output: Each step has its own condor job template file (.jtf) and output see below:

Steps:
    Trimmomatic         -- http://www.usadellab.org/cms/?page=trimmomatic
    Fastqc              -- http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
    bowtie2             -- http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
    bwa mem             -- http://bio-bwa.sourceforge.net/
    bam2wig             -- http://search.cpan.org/~tjparnell/Bio-ToolBox-1.24001/lib/Bio/ToolBox.pm
    picard/CleanSam.jar -- http://broadinstitute.github.io/picard/
    samtools            -- http://www.htslib.org/
    HTSeq               -- http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html
    Normalization(RPKM) -- /home/GLBRCORG/mplace/scripts/RPKM.py 

Steps:
*******************************************************************************
Trimmomatic:
        HEADCROP       = 5  
        LEADING        = 3 
        TRAILING       = 3
        SLIDINGWINDOW  = 3:30   (window = 3, min avg quality for window = 30)
        MINLEN         = 36        
        -phred33 
        -threads 8
        SE is single-end
        PE is paired-end
        -trimlog <logfile>

    Single-end
    java -jar <path to trimmomatic jar> SE [-threads <threads>] -phred33 
    [-trimlog <logFile>] <input> <output> 

    on GLBRC scarcity:
    java -jar  /opt/bifxapps/Trimmomatic-0.30/trimmomatic-0.30.jar

    OUTPUT: trimmed fastq files
    
*******************************************************************************
FastQC:
    To run non-interactively simply a list of files to process on the commandline

    fastqc somefile.txt someotherfile.txt

    You can specify as many files to process in a single run as you like

    If you want to save your reports in a folder other than the folder 
    
    --outdir=/some/other/dir/
    -quiet   = only report errors
    
    on GLBRC scarcity:
    /opt/bifxapps/FastQC/fastqc --help
    
    OUTPUT: QC report after trimming with Trimmomatic

*******************************************************************************
MAPPING:  default is Bowtie2, but user may choose bwa mem
*******************************************************************************
Bowtie2:
    -p 8           # specified number of parallel search threads
    --phred33 
    -N 1           # Sets the number of mismatches to allowed in a seed alignment
                     during multiseed alignment. Can be set to 0 or 1.
    -x $REF        # index reference

    -U $READS      # read
    -S $OUT.sam    # output sam file

    on GLBRC scarcity:
    /opt/bifxapps/bin/bowtie2
    
    OUTPUT: sam file

Bwa mem:
    -t 8          # number of threads
    -M $REFERENCE # reference genome file
  example:
    bwa mem -t 8 -M $REFERENCE $file 1> $out.sam 2>>$OUTPUT_DIR/Bwa_run.log
    
    OUTPUT: sam file

*******************************************************************************
picard

 Clean the SAM file
    This soft-clips an alignment that hangs off the end of its reference sequence.
    This will print out all the errors that it ignores (MAPQ errors)

java -Xmx15g -jar /opt/bifxapps/picard/CleanSam.jar I=$OUT.sam O=$OUT.cleaned.sam
rm $OUT.sam

 Add the RG header and sort the SAM file
    This will print out all the errors that it ignores (MAPQ errors)

java -Xmx15g -jar /opt/bifxapps/picard/AddOrReplaceReadGroups.jar I=$OUT.cleaned.sam 
O=$OUT.final.sam SO=coordinate LB=$REF.fasta PL=ILLUMINA PU=unknown SM=$OUT 
VALIDATION_STRINGENCY=LENIENT

rm $OUT.cleaned.sam
    on GLBRC scarcity:
    java -Xmx15g -jar /opt/bifxapps/picard/AddOrReplaceReadGroups.jar
    
    OUTPUT: sam file
*******************************************************************************
samtools

Make the BAM file, sort and index it
    samtools view -uS -t $REF.fasta.fai $OUT.final.sam | samtools sort - $OUT.sorted

samtools index $OUT.sorted.bam 

    on GLBRC scarcity:
    /opt/bifxapps/bin/samtools
    
    OUTPUT: sorted bam file and index for bam file
*******************************************************************************
bam2wig.pl

This script will convert alignments from a Bam file into enumerated point data
in a wig format.

bam2wig.pl --in bamFile --pos mid --strand --rpm --out

This works on  scarcity-1, scarcity-5, scarcity-6 but may fail elsewhere

    OUTPUT:  gzipped wig file

*******************************************************************************
HTSeq

Given a file with aligned sequencing reads and a list of genomic features, 
count how many reads map to each feature.

htseq-count -t CDS -i Parent samFile  gff 

    OUTPUT: htseq text file

*******************************************************************************
RPKM

RPKM normalization - /home/GLBRCORG/mplace/scripts/RPKM.py

RPKM.py -d <directory> -g REF.gff [-f gene]

Run RPKM on all HTSeq files in Current Directory.
Requires bam files to be in the same directory.
Header of bam files are used to get genome size.

    OUTPUT: RPKM text file , one column per sample
    
*******************************************************************************
@Required:
    Picard tools
    bowtie2
    fastqc
    trimmomatic
    python (snakemake)
    bowtie2-build <reference.in> <basename for index files>
    samtools faidx <referenc.in>

@Author: mplace
   
Dependencies: Python 3
              pydagman modules provided by Branden Timm

Date:   1/05/2016

pydagman located: /home/GLBRCORG/mplace/anaconda3/lib/python3.4/site-packages/pydagman
"""
import argparse
import os
import re
import socket
import sys
from pydagman.dagfile import Dagfile
from pydagman.job import Job

ref = { 'R64' : ( "/home/GLBRCORG/mplace/data/reference/S288C_reference_genome_R64-1-1_20110203/s.cerevisiae-R64-1-1",
                  "/home/GLBRCORG/mplace/data/reference/S288C_reference_genome_R64-1-1_20110203/S288C_reference_sequence_R64-1-1_20110203.fasta",
                  "/home/GLBRCORG/mplace/data/reference/S288C_reference_genome_R64-1-1_20110203/saccharomyces_cerevisiae_R64-1-1_20110208_noFasta.gff",
                  "/home/GLBRCORG/mplace/data/reference/S288C_reference_genome_R64-1-1_20110203/S288C_reference_sequence_R64-1-1_20110203.fsa.fai",
                  "/home/GLBRCORG/mplace/data/reference/S288C_reference_genome_R64-1-1_20110203/S288C_reference_sequence_R64-1-1_20110203.fasta"),
        
        'Y22' : ( "/home/GLBRCORG/mplace/data/reference/Y22-3/Y22-3-bowtie",
                  "/home/GLBRCORG/mplace/data/reference/Y22-3/Y22-3.fasta",
                  "/home/GLBRCORG/mplace/data/reference/Y22-3/Y22-3_Final_GFF.gff",
                  "/home/GLBRCORG/mplace/data/reference/Y22-3/Y22-3.fasta.fai",
                  "/home/GLBRCORG/mplace/data/reference/Y22-3/Y22-3.fasta" ),
        
        'PAN' : ( "/home/GLBRCORG/mplace/data/reference/PanGenome/PanGenome-Final-R64-2-1-bowtie",
                  "/home/GLBRCORG/mplace/data/reference/PanGenome/PanGenome-Final-R64-2-1.fasta",
                  "/home/GLBRCORG/mplace/data/reference/PanGenome/PanGenome-Final-R64-2-1.gff",
                  "/home/GLBRCORG/mplace/data/reference/PanGenome/PanGenome-Final-R64-2-1.fasta.fai",
                  "/home/GLBRCORG/mplace/data/reference/PanGenome/PanGenome-Final-R64-2-1.fasta" ) 
    }

def trimCondorFile():
    """
    Step 1
    Create condor job template file to run trimmomatic.

    Command:
    java -Xmx6g -jar ~/bin/trimmomatic SE -phred33 -trimlog trimlog.out test.fastq.gz
    trimmed.out.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:3:30 MINLEN:36

    To use with condor you must know the main class.
    extract jar archive: jar xf trimmomatic-0.30.jar

    look in the META-INF/MANIFEST.MF for main class

    Manifest-Version: 1.0
    Ant-Version: Apache Ant 1.8.3
    Created-By: 1.7.0_09-icedtea-mockbuild_2013_01_14_23_04-b00 (Oracle Corporation)
    Main-Class: org.usadellab.trimmomatic.Trimmomatic

    pass the main class in the Arguments line.    
  
    """
    with open('trimCondor.jtf', 'w') as submit:
        submit.write( "Universe                 = java\n" )
        submit.write( "Executable               = /opt/bifxapps/Trimmomatic-0.30/trimmomatic-0.30.jar\n" )
        submit.write( "jar_files                = /opt/bifxapps/Trimmomatic-0.30/trimmomatic-0.30.jar\n" )
        submit.write( "java_vm_args             = -Xmx6g\n" )
        submit.write( "Arguments                = org.usadellab.trimmomatic.Trimmomatic SE -phred33 $(fastq) $(outfile) LEADING:3 TRAILING:3 SLIDINGWINDOW:3:30 MINLEN:36\n" )
        submit.write( "Notification             = Never\n" )
        submit.write( "Should_Transfer_Files    = Yes\n" )
        submit.write( "When_To_Transfer_Output  = On_Exit\n" )
        submit.write( "Transfer_Input_Files     = $(fastq)\n" )
        submit.write( "Error                    = $(job).submit.err\n" )
        submit.write( "Log                      = $(job).submit.log\n" )
        submit.write( "request_memory           = 8000\n" )
        submit.write( "request_disk             = 400000\n" )
        submit.write( "Queue\n" )
    submit.close()        

def fastqcCondorFile():
    """
    NOT USED YET, need to find a way to call the executable, command
    below does not work.
    Step 2
    Create condor job template file to run FastQC.
    This provides some quality control checks on trimmed sequence data.

    Command:
    /opt/bifxapps/FastQC/fastqc input.fastq
    """
    with open('fastqcCondor.jtf', 'w') as submit: 
        submit.write( "Universe                 = java\n" )
        submit.write( "Executable               = /opt/bifxapps/FastQC/fastqc\n" )
        submit.write( "Arguments                = $(read)\n" )
        submit.write( "Notification             = Never\n" )
        submit.write( "Should_Transfer_Files    = Yes\n" )
        submit.write( "When_To_Transfer_Output  = On_Exit\n" )
        submit.write( "Transfer_Input_Files     = $(read)\n" )
        submit.write( "Error                    = $(job).submit.err\n" )
        submit.write( "Log                      = $(job).submit.log\n" )
        submit.write( "request_memory           = 20000\n" )
        submit.write( "request_disk             = 5000000\n" )
        submit.write( "Queue" )
    submit.close()

def bowtie2CondorFile():
    """
    Create condor job template file to run bowtie2.
    """
    pass


def bwaCondorFile():
    """
    Create condor job template file to run bwa mem.
    """
    with open('bwaCondor.jtf', 'w') as submit: 
        submit.write( "Universe                 = vanilla\n" )
        submit.write( "Executable               = /opt/bifxapps/bin/bwa\n" )
        submit.write( "Arguments                = mem -t 8 -M $(reference) $(read)\n" )
        submit.write( "Notification             = Never\n" )
        submit.write( "Should_Transfer_Files    = Yes\n" )
        submit.write( "When_To_Transfer_Output  = On_Exit\n" )
        submit.write( "Transfer_Input_Files     = $(reference),$(read)\n" )
        submit.write( "Output			 = $(outfile)\n" )
        submit.write( "Error                    = $(job).submit.err\n" )
        submit.write( "Log                      = $(job).submit.log\n" )
        submit.write( "request_memory           = 20000\n" )
        submit.write( "request_disk             = 5000000\n" )
        submit.write( "Queue" )
    submit.close()

def cleanSamFile():
    """
    Create condor job template file to run picard tools CleanSam on bwa results.
    """
    with open('cleanSam.jtf', 'w') as submit:
        submit.write( "Universe                 = java\n" )
        submit.write( "Executable               = /opt/bifxapps/picard/picard.jar\n" )
        submit.write( "jar_files                = /opt/bifxapps/picard/picard.jar\n" )
        submit.write( "java_vm_args             = -Xmx8g\n" )
        submit.write( "Arguments                = picard.cmdline.PicardCommandLine CleanSam I=$(sam) O=$(outfile) QUIET=true VERBOSITY=null TMP_DIR=/tmp\n" )
        submit.write( "Notification             = Never\n" )
        submit.write( "Should_Transfer_Files    = Yes\n" )
        submit.write( "When_To_Transfer_Output  = On_Exit\n" )
        submit.write( "Transfer_Input_Files     = $(sam)\n" )
        submit.write( "Error                    = $(job).submit.err\n" )
        submit.write( "Log                      = $(job).submit.log\n" )
        submit.write( "request_memory           = 8000\n" )
        submit.write( "request_disk             = 400000\n" )
        submit.write( "Queue" )
    submit.close()

def addReadGpSam(ref):
    """
    Create condor job template file to add RG header to
    sam file using Picard tools.
    """
    with open('addReadGpSam.jtf', 'w') as submit:
        submit.write( "Universe                 = java\n" )
        submit.write( "Executable               = /opt/bifxapps/picard/picard.jar\n" )
        submit.write( "jar_files                = /opt/bifxapps/picard/picard.jar\n" )
        submit.write( "java_vm_args             = -Xmx8g\n" )
        submit.write( "Arguments                = picard.cmdline.PicardCommandLine AddOrReplaceReadGroups I=$(cleanSam) O=$(outfile) SO=coordinate LB=" + ref + " PL=ILLUMINA PU=unknown SM=$(fastq) VALIDATION_STRINGENCY=LENIENT TMP_DIR=/tmp\n" )
        submit.write( "Notification             = Never\n" )
        submit.write( "Should_Transfer_Files    = Yes\n" )
        submit.write( "When_To_Transfer_Output  = On_Exit\n" )
        submit.write( "Transfer_Input_Files     = $(cleanSam)\n" )
        submit.write( "Error                    = $(job).submit.err\n" )
        submit.write( "Log                      = $(job).submit.log\n" )
        submit.write( "request_memory           = 8000\n" )
        submit.write( "request_disk             = 400000\n" )
        submit.write( "Queue" )
    submit.close()  

def samToBamFile():
    """
    Convert sam file, output of addReadGpSam(), to bam for sorting by samtools.
    """
    with open('samToBam.jtf', 'w') as submit:
        submit.write( "Universe                 = vanilla\n" )
        submit.write( "Executable               = /opt/bifxapps/bin/samtools\n" )
        submit.write( "Arguments                = view -bS -t $(reference) -o $(bam) $(sam)\n" )
        submit.write( "Notification             = Never\n" )
        submit.write( "Should_Transfer_Files    = Yes\n" )
        submit.write( "When_To_Transfer_Output  = On_Exit\n" )
        submit.write( "Transfer_Input_Files     = $(reference),$(sam)\n" )
        submit.write( "Error                    = $(job).submit.err\n" )
        submit.write( "Log                      = $(job).submit.log\n" )
        submit.write( "request_memory           = 20000\n" )
        submit.write( "request_disk             = 5000000\n" )
        submit.write( "Queue" )
    submit.close()
        
def sortSamFile():
    """
    Use samtools to Sort bam file, takes output from samToBamFile.
    """
    with open('sortSam.jtf', 'w') as submit:
        submit.write( "Universe                 = vanilla\n" )
        submit.write( "Executable               = /opt/bifxapps/bin/samtools\n" )
        submit.write( "Arguments                = sort $(bam) $(out)\n" )
        submit.write( "Notification             = Never\n" )
        submit.write( "Should_Transfer_Files    = Yes\n" )
        submit.write( "When_To_Transfer_Output  = On_Exit\n" )
        submit.write( "Transfer_Input_Files     = $(bam)\n" )
        submit.write( "Error                    = $(job).submit.err\n" )
        submit.write( "Log                      = $(job).submit.log\n" )
        submit.write( "request_memory           = 20000\n" )
        submit.write( "request_disk             = 5000000\n" )
        submit.write( "Queue" )
    submit.close()

def indexBamFile():
    """
    Use samtools to index bam file, takes output from sortSamFile.
    """
    with open('indexBam.jtf', 'w') as submit:
        submit.write( "Universe                 = vanilla\n" )
        submit.write( "Executable               = /opt/bifxapps/bin/samtools\n" )
        submit.write( "Arguments                = index $(bam)\n" )
        submit.write( "Notification             = Never\n" )
        submit.write( "Should_Transfer_Files    = Yes\n" )
        submit.write( "When_To_Transfer_Output  = On_Exit\n" )
        submit.write( "Transfer_Input_Files     = $(bam)\n" )
        submit.write( "Error                    = $(job).submit.err\n" )
        submit.write( "Log                      = $(job).submit.log\n" )
        submit.write( "request_memory           = 20000\n" )
        submit.write( "request_disk             = 5000000\n" )
        submit.write( "Queue" )
    submit.close()

def htSeqFile( strandedness ):
    """
    Call HTSeq-count to get gene counts for each sample.
    -s parameter indicates reverse strand library construction
    """
    with open('htseq.jtf', 'w') as submit:
        submit.write( "Universe                 = vanilla\n" )
        submit.write( "Executable               = /opt/bifxapps/python/bin/htseq-count\n" )
        #submit.write( "Environment              = \"PYTHONPATH=/opt/bifxapps/python/lib64/python2.6/site-packages/:/opt/bifxapps/python/lib/python2.6/site-packages/:/home/GLBRCORG/mplace/anaconda3/lib/python3.4/site-packages\"")
        submit.write( "getenv                   = True\n")
        if strandedness == 1:
            submit.write( "Arguments                = -f bam -t CDS -i Parent -s reverse $(bam) $(gff)\n" )
        else:
            submit.write( "Arguments                = -f bam -t CDS -i Parent $(bam) $(gff)\n" )
        submit.write( "Notification             = Never\n" )
        submit.write( "Should_Transfer_Files    = Yes\n" )
        submit.write( "When_To_Transfer_Output  = On_Exit\n" )
        submit.write( "Transfer_Input_Files     = $(bam)\n" )
        submit.write( "output                   = $(out)\n" )
        submit.write( "Error                    = $(job).submit.err\n" )
        submit.write( "Log                      = $(job).submit.log\n" )
        submit.write( "request_memory           = 20000\n" )
        submit.write( "request_disk             = 5000000\n" )
        submit.write( "Queue" )
    submit.close()
    
def rpkmFile():
    """
    Run RPKM normalization on all HTSeq files in current working directory.
    
    Calls /home/GLBRCORG/mplace/scripts/RPKM.py

    RPKM.py -d <directory> -g REF.gff

    The assumption is that as part of this pipeline HTSeq has been run
    using CDS as the genome feature.
    """
    with open('rpkm.jtf', 'w') as submit:
        submit.write( "Universe                 = vanilla\n" )
        submit.write( "Executable               = /home/GLBRCORG/mplace/scripts/RPKM.py\n" )
        submit.write( "Arguments                = -d $(cwd) -f $(gff)\n" )
        submit.write( "Notification             = Never\n" )
        submit.write( "Should_Transfer_Files    = Yes\n" )
        submit.write( "When_To_Transfer_Output  = On_Exit\n" )
        submit.write( "Transfer_Input_Files     = $(gff)\n" )
        submit.write( "Error                    = $(job).submit.err\n" )
        submit.write( "Log                      = $(job).submit.log\n" )
        submit.write( "request_memory           = 20000\n" )
        submit.write( "request_disk             = 5000000\n" )
        submit.write( "Queue" )
    submit.close()

def bamToWigFile():
    """
    Create wig file for viewing alignments in mochi view
    """
    pass

def main():
    """
    main() 
    """
    cmdparser = argparse.ArgumentParser( description="RNA-Seq Pipeline, condor version",
                                         usage='%(prog)s -f <fastq file list.txt> [optional args: -a -r -d -ref ]', prog='rnaSeqCondor.py')
    cmdparser.add_argument('-a', '--aligner', action='store',      dest='ALIGNER', help='Default aligner is Bowtie2, to use Bwa mem: -a bwamem')
    cmdparser.add_argument('-d', '--detail',  action='store_true', dest='DETAIL',  help='Print a more detailed description of program.')    
    cmdparser.add_argument('-r', '--reverse', action='store_true', dest='REVERSE', help='HTSeq -s reverse, for Biotech GEC data, optional.')
    cmdparser.add_argument('-ref', '--reference', action='store',  dest='REFERENCE', help='Choose reference -ref R64 (SGD R64-1-1) -ref Y22-3 (GLBRC)' )
    cmdResults = vars(cmdparser.parse_args())

    fastq = []              # List of fastq files to process
    cwd   = os.getcwd()     # Current working directory

    #if no args print help
    if len(sys.argv) == 1:
        print("")
        cmdparser.print_help()
        sys.exit(1)

    if cmdResults['DETAIL']:
        print('\nrnaSeqCondor.py')
        print("\nPurpose: Condor implementation of the Gasch lab RNA-Seq pipeline.")
        print("\t Maybe be run using GLOW or from command line.")
        print("\nInput : none, just run in directory with fastq files to process.")
        print("\tGlow will handle the job submission automatically.\n")
        print("This script must be run on scarcity's condor submit node, scarcity-cm.glbrc.org " )
        print("\nPlease use a dedicated directory for running pipeline.")
        print("If your are running this manually do the following:\n")
        print("\tCreate your directory and copy or link your fastq files into that directory.\n")
        print("To run default enter:  /home/GLBRCORG/mplace/scripts/rnaSeqCondor.py \n")
        print("Optional Parameters:") 
        print("\t-r  this will use \"-s reverse\" parameter for HTSeq.\n")
        print("\t-a  bwamem changes aligner from default bowtie2 to bwa mem\n")
        print("\t-ref  change default reference, usage:  -ref Y22-3")
        print("\t    Current Reference List:")
        print("\t\tR64-1-1 -- default equivalent to UCSC sacCer3" )
        print("\t\tY22-3   -- GLBRC assembly")
        print("\nOutline of steps & commands used in pipeline:\n")
        print("  1) Trimmomatic ")
        print("\t/opt/bifxapps/Trimmomatic-0.30/trimmomatic-0.30.jar SE -phred33 input_fastq ")
        print("\toutFile LEADING:3 TRAILING:3 SLIDINGWINDOW:3:30 MINLEN:36\n")
        print("  2) Fastqc ")
        print("\t/opt/bifxapps/FastQC/fastqc trimmed_fastq\n")
        print("  3) Alignment ")
        print("\tBowtie2 -p 8 --phred33 -N 1 -x  referenceFile -U  trimmed_fastq -S outFile" )
        print("\tor")
        print("\tbwa mem -t 8 -M referenceFile trimmed_fastq\n")
        print("  4) Picard-tools ")
        print("\t/opt/bifxapps/picard-tools/CleanSam.jar I= samFile O=outFile")
        print("\t/opt/bifxapps/picard-tools/AddOrReplaceReadGroups.jar 'I=', inFile ")
        print("\tO=outFile SO=coordinate LB=S288C_reference_sequence_R64-1-1_20110203.fasta\n\tPL=ILLUMINA PU=unk RGSM=SampleName\n")
        print("  5) samtools ")
        print("\tsam to bam: samtools view -bS -t reference.fsa.fai -o bam samFile")
        print("\t  sort bam: samtools sort input_bam sorted_Bam ")
        print("\t index bam: samtools index sorted_Bam\n")
        print("  6) HTSeq ")
        print("\t/opt/bifxapps/python/bin/htseq-count -t CDS -i Parent inputFile Ref_gff")
        print("\t -- if '-r' specfied -s reverse will be used as well\n")
        print("  7) RPKM ")
        print("\t/home/GLBRCORG/mplace/bin/Normalize.jar dir Ref_gff RPKM.results --gene\n")
        print("  8) bam2wig.pl ")
        print("\t/opt/bifxapps/biotoolbox/scripts/bam2wig.pl --in bamFile --pos mid --strand --rpm --out outFile\n")
        print("  9) findreplace_WIG.pl")
        print("\t/home/GLBRCORG/mplace/scripts/findreplace_WIG.pl\n")
        print("The results are organized into the following subdirectories:")
        print("\t alignments/ fastq/ fastqc/ htseq/ log/ wig/ ")
        print("\t RPKM results are written to a file called: RPKM.results\n")
        print("This script is designed to run on the GLBRC condor submit node scarcity-cm")
        print("See Mike Place for problems with this script.")
        sys.exit(1)

    hostname = socket.gethostname()
    if hostname != 'scarcity-cm.glbrc.org':
        print("Program must be run on scarcity's condor submit node.")
        print("\tssh to scarcity-cm.glbrc.org")
        sys.exit(1)

    # Get reference genome to use
    if cmdResults['REFERENCE'] == 'Y22-3' or cmdResults['REFERENCE'] == 'y22-3':
        reference = 'Y22'
    else:
        reference = 'R64'

    # Get aligner to use
    if cmdResults['ALIGNER'] is not None:
        if cmdResults['ALIGNER'] == 'bwamem'
        aligner = 'bwamem'
    else:
        aligner='bowtie2'

    # Get input file listing fastq files to process
    for f in os.listdir(cwd):
        if re.search(r'.fastq',f):
            fastq.append(f)

    with open('bwamem.log','w') as log:
        log.write("Running bwa mem on condor\n")
        log.write("Using the following input files:\n")
        for file in fastq:
            log.write('\t%s\n' %(file))
        log.write("Reference : %s\n" %(reference))
        log.write("\n")
        log.close()
        
    cwd = os.getcwd()
    
    #create Dagfile object
    mydag = Dagfile()
    numJobs = len(fastq)

    num = 1
    for f in fastq:
        
        trimJob = Job('trimCondor.jtf', 'job' + str(num))   # set trimmomatic job file
        mydag.add_job(trimJob)
        trimJob.pre_skip("1")
        trimJob.add_var('job', 'job' + str(num))             # setup variable to substitue
        trimJob.add_var('fastq', f )
        trimName = re.sub(r"fastq","trim.fastq", f)
        trimJob.add_var('outfile', trimName)  
        num += 1   
        
        bwaJob =  Job('bwaCondor.jtf', 'job' + str(num))     # set bwa job file
        mydag.add_job(bwaJob)   
        bwaJob.pre_skip("1")
        bwaJob.add_var('job', 'job' + str(num))              # setup variable to substitue
        bwaJob.add_var('reference', ref[reference][1])
        bwaJob.add_var('read', f )
        outName = re.sub(r"fastq","sam", trimName)
        bwaJob.add_var('outfile', outName)
        bwaJob.add_parent(trimJob)
        num += 1
          
        cleanJob =    Job('cleanSam.jtf', 'job' + str(num))  # set up clean sam job file
        cleanJob.pre_skip("1")
        cleanJob.add_var('job', 'job' + str(num))
        cleanJob.add_var( 'sam', outName )
        cleanName = re.sub( r"fastq", "clean.sam", f)
        cleanJob.add_var( 'outfile', cleanName)
        parent = 'job' + str(num)
        cleanJob.add_parent(bwaJob)                           # Make bwa job parent of clean sam job
        mydag.add_job(cleanJob)
        num += 1

        readGpJob = Job('addReadGpSam.jtf', 'job' + str(num)) # set up AddOrReplaceReadGroups job
        readGpJob.pre_skip("1")
        readGpJob.add_var('job', 'job' + str(num))
        readGpJob.add_var('cleanSam', cleanName)
        readGpJob.add_var('fastq', f)
        outfile = re.sub(r"fastq", "final.sam", f)
        readGpJob.add_var('outfile', outfile)
        parent = 'job' + str(num)
        readGpJob.add_parent(cleanJob)
        mydag.add_job(readGpJob)
        num += 1

        samToBamJob  = Job('samToBam.jtf', 'job' + str(num))   # set up sam to bam job
        samToBamJob.pre_skip("1")
        samToBamJob.add_var('job', 'job' + str(num))
        samToBamName = re.sub(r"fastq", "final.sam", f)
        samToBamJob.add_var('sam', samToBamName)
        bamName      = re.sub(r"sam", "bam", samToBamName)
        samToBamJob.add_var('bam', bamName)
        samToBamJob.add_var('reference', ref[reference][3])
        parent = 'job' + str(num)
        samToBamJob.add_parent(readGpJob)
        mydag.add_job(samToBamJob)
        num += 1

        sortSamJob = Job('sortSam.jtf', 'job' + str(num))      # set up sort sam job
        sortSamJob.pre_skip("1")
        sortSamJob.add_var('job', 'job' + str(num))
        sortSamJob.add_var('bam', bamName )
        sortName = re.sub(r"bam", "sort", bamName)
        sortSamJob.add_var('out', sortName )
        parent = 'job' + str(num)
        sortSamJob.add_parent(samToBamJob)
        mydag.add_job(sortSamJob)
        num += 1

        indexBamJob = Job('indexBam.jtf', 'job' + str(num))     # set up index sam job
        indexBamJob.pre_skip("1")
        indexBamJob.add_var('job', 'job' + str(num))
        indexBamJob.add_var('bam', sortName + '.bam' )
        parent = 'job' + str(num)
        indexBamJob.add_parent(sortSamJob)
        mydag.add_job(indexBamJob)
        num += 1

        htseqJob = Job('htseq.jtf', 'job' + str(num))           # set up HTSeq job
        htseqJob.pre_skip("1")
        htseqJob.add_var('job', 'job' + str(num))
        htseqJob.add_var('bam', sortName + '.bam' )
        htseqJob.add_var('gff', ref[reference][2] )             
        htseqName = re.sub('.sam', '_HTseqOutput.txt', samToBamName)
        htseqJob.add_var('out', htseqName )
        parent = 'job' + str(num)
        htseqJob.add_parent(indexBamJob)
        mydag.add_job(htseqJob)
        num += 1

        rpkmJob = Job('rpkm.jtf', 'job' + str(num))             # set up RPKM job
        rpkmJob.pre_skip("1")
        rpkmJob.add_var('job', 'job' + str(num))
        rpkmJob.add_var('cwd', cwd)
        rpkmJob.add_var('gff', ref[reference][2])
        parent = 'job' + str(num)
        mydag.add_job(rpkmJob)
        num += 1

        
    # write trimmomatic submit file
    trimCondorFile()
    # write fastqc submit file
    fastqcCondorFile()
    # write bwa submit file
    bwaCondorFile()
    # write clean sam submit file
    cleanSamFile()
    # write add read group to sam condor 
    addReadGpSam(reference)
    # write sam To Bam submit file
    samToBamFile()
    # write sort sam submit file
    sortSamFile()
    # write index bam submit file
    indexBamFile()
    # write HTSeq submit file
    if cmdResults['REVERSE'] == 1:
        strandedness = 1
    else:
        strandedness = 0
    htSeqFile(strandedness)

    
    
    mydag.save('MasterDagman.dsf')
     
    # clean up unnessary files
    # remove original sam file, generated w/ bwa
    
      

if __name__ == "__main__":
    main()

