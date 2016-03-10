# reference.py provides single point of access for references uses with
# the HTcondor RNA-Seq pipeline.

# reference dictionary to stores locations of reference fasta, gff, dict files
# key = reference name, value is a tuple where the order is defined as:
#   [0] = bowtie2 reference
#   [1] = bwa mem reference
#   [2] = gff file
#   [3] = samtools index
#   [4] = picard used for nameing
# default reference = R64 (SGD R64-1-1 = UCSC sacCer3)
# Y22 = reference of S. cerevisiae Y22-3 GLBRC sequenced strain
# PAN = Pan-genome, composed of R64-2-1, sequence regions from Strope et. al
#       , novel genes from Borneman strains and JAY291 (Argueso et.al)

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
