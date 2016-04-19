rnaSeqCondor.py
  
Purpose:

    HTCondor implementation of the Gasch lab RNA-Seq pipeline.
    Pipeline can be run 2 ways, from the command line or from
    the GLOW website.

    Fastq files should be gzipped and listed in put files with full path.

    Job implemented as a dagman using Branden Timm's Python Pydagman modules.


To get a detailed description run: 

    ./rnaSeqCondor.py -d


Outline of steps & commands used in pipeline:

  1) Trimmomatic version 0.3 

          /opt/bifxapps/Trimmomatic-0.30/trimmomatic-0.30.jar SE -phred33 input_fastq 
          outFile LEADING:3 TRAILING:3 SLIDINGWINDOW:3:30 MINLEN:36

  2) Fastqc version 0.10.1

          /opt/bifxapps/FastQC/fastqc trimmed_fastq

  3) Alignment 

          Bowtie2 version 2.2.2

          Bowtie2 -p 8 --phred33 -N 1 -x  referenceFile -U  trimmed_fastq -S outFile
                                              
          or

          bwa version 0.7.12-r1039

          bwa mem -t 8 -M referenceFile trimmed_fastq

  4) Picard-tools version 1.98(1547)

          /opt/bifxapps/picard-tools/CleanSam.jar I= samFile O=outFile

          /opt/bifxapps/picard-tools/AddOrReplaceReadGroups.jar 'I=', inFile 

          O=outFile SO=coordinate LB=S288C_reference_sequence_R64-1-1_20110203.fasta

          PL=ILLUMINA PU=unk RGSM=SampleName

  5) samtools version 1.2 (using htslib 1.2.1)

          sam to bam: samtools view -bS -t reference.fsa.fai -o bam samFile

          sort bam: samtools sort input_bam sorted_Bam 

          index bam: samtools index sorted_Bam

  6) HTSeq version 0.6.0

         /opt/bifxapps/python/bin/htseq-count -t CDS -i Parent inputFile Ref_gff

         -- if '-r' specfied -s reverse will be used as well

  7) RPKM version 1.0

         /home/GLBRCORG/mplace/scripts/RPKM.py
                                                                                                                                                      
         usage: RPKM.py -r refname [-f gene] 

         Run RPKM on all HTSeq files in Current Directory.

  8) bam2wig.pl version 1.12.5

         /opt/bifxapps/biotoolbox/scripts/bam2wig.pl --in bamFile --pos mid --strand --rpm --out outFile


  The results are organized into the following subdirectories:

        fastqc/         -- Fastqc results, zipped
                                                                             
        htseq/          -- HTSeq results 
                                                                                             
        log/            -- log and error files
                                                                                                             
        wig/            -- wig files for visualization
                                                                                                                             
        RPKM results are written to a file called: RPKM.results

        All alignment files (.sam & .bam) are removed to save disk space.


Parameters:

    -a  default aligner bowtie2, to use bwa : -a  bwamem

    -f  input file listing full path to fastq.gz files one per line

    -i  GLOW WorkFlow ID number (required).

    -r  this will use "-s reverse" parameter for HTSeq.

    -ref  change default reference, usage:  -ref Y22-3
    
        Current Reference List:
        R64-1-1 -- default equivalent to UCSC sacCer3
        R64-2-1 -- Most recent SGD S.cerevisiae genome reference
        PAN     -- S.cerevisiae PanGenome reference 
        Y22-3   -- GLBRC strain Y22-3 S.cerevisiae reference

*******************************************************************************
RPKM.py

    Calculate RPKM for all HTSeq results in directory.  Bam files are 
    expected to be present as well.

    usage: RPKM.py -r genome [-f feature]

    example: ./RPKM.py  -r R64 

    Run RPKM on all HTSeq files in Current Directory.

    optional arguments:
    -h, --help       show this help message and exit
    -f, --feature   Feature type to use w/ GFF file, default = CDS.
    -r, --reference   Genome Name [R64, R64-1, PAN, Y22]
    -i, --info       Detailed description of program.
                            
******************************************************************************
rnaSeqCleanUp.py

    Purpose: Post Condor RNA-Seq pipeline run after all other condor jobs are
             complete in current directory.
                    
    Input:   Genome [R64, R64-1, PAN, Y22] 
                              
                                                   
    Output: The following directories are created fastqc, htseq,log, wig 
            and the associated files are moved into them. 
            RPKM is run on all HTSeq files with results in RPKM.results file.

******************************************************************************

