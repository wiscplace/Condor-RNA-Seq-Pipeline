rnaSeqCondor.py
  
Purpose:

    HTCondor implementation of the Gasch lab RNA-Seq pipeline.

    All fastq files in the working directory will be processed.

    Job implemented as a dagman using Branden Timm's Python Pydagman modules.

Outline of steps & commands used in pipeline:

  1) Trimmomatic 

          /opt/bifxapps/Trimmomatic-0.30/trimmomatic-0.30.jar SE -phred33 input_fastq 
          outFile LEADING:3 TRAILING:3 SLIDINGWINDOW:3:30 MINLEN:36

  2) Fastqc 

          /opt/bifxapps/FastQC/fastqc trimmed_fastq

  3) Alignment 

          Bowtie2 -p 8 --phred33 -N 1 -x  referenceFile -U  trimmed_fastq -S outFile
                                              
          or

          bwa mem -t 8 -M referenceFile trimmed_fastq

  4) Picard-tools 

          /opt/bifxapps/picard-tools/CleanSam.jar I= samFile O=outFile

          /opt/bifxapps/picard-tools/AddOrReplaceReadGroups.jar 'I=', inFile 

          O=outFile SO=coordinate LB=S288C_reference_sequence_R64-1-1_20110203.fasta

          PL=ILLUMINA PU=unk RGSM=SampleName

  5) samtools 

          sam to bam: samtools view -bS -t reference.fsa.fai -o bam samFile

          sort bam: samtools sort input_bam sorted_Bam 

          index bam: samtools index sorted_Bam

  6) HTSeq 

         /opt/bifxapps/python/bin/htseq-count -t CDS -i Parent inputFile Ref_gff

         -- if '-r' specfied -s reverse will be used as well

  7) RPKM 

         /home/GLBRCORG/mplace/scripts/RPKM.py
                                                                                                                                                      
         usage: RPKM.py -r refname [-f gene] 

         Run RPKM on all HTSeq files in Current Directory.

  8) bam2wig.pl 

         /opt/bifxapps/biotoolbox/scripts/bam2wig.pl --in bamFile --pos mid --strand --rpm --out outFile


  The results are organized into the following subdirectories:

        fastqc/         -- Fastqc results, zipped
                                                                             
        htseq/          -- HTSeq results 
                                                                                             
        log/            -- log and error files
                                                                                                             
        wig/            -- wig files for visualization
                                                                                                                             
        RPKM results are written to a file called: RPKM.results

        All alignment files (.sam & .bam) are removed to save disk space.

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

