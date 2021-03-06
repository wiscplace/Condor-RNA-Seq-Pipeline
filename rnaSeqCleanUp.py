#!/home/GLBRCORG/mplace/anaconda3/bin/python
"""
@Program: rnaSeqCleanUp.py

@Purpose: Post Condor RNA-Seq pipeline DAGman final process and clean up
          
@Input:  Reference genome to use, options are R64, R64-2, PAN, Y22            
                     
@Output: The following directories are created fastqc, htseq, wig, log
         and associate files are moved there. RPKM is calculated for all samples
         output is in RPKM.results.

@author: Mike Place
@Date:   3/1/2016

TODO:
get submitter name
use submitter name to create dir on bigdata/processed_data

"""
import os
import re
import shutil
import subprocess
import sys
import reference as r
import sendEmail as mail

def cleanUp( cwd, submitter, rnaDir, workflowID, token):
    """
    Delete unneeded files and move output files to appropriate directory on big data
    files types deleted: bam, fastq, sam
    directories created for log, wig, htseq files, these files remain in users home directory
    A copy of the RPKM.results file will be put in the users directory on big data.
    GLOW will be updated with the results and job completion.
    """
    # copy RPKM.results file to bigdata
    os.mkdir('/mnt/bigdata/processed_data/' + submitter + '/' + rnaDir)
    shutil.copy('RPKM.results', '/mnt/bigdata/processed_data/' + submitter + '/' + rnaDir)
    with open('pipeline.log', 'a') as log:
        log.write('\n RPKM.results file written to:\n')
        log.write('/mnt/bigdata/processed_data/' + submitter + '/' + rnaDir)
    log.close()
    
    cwd = cwd + "/"
    # remove bam files
    [ os.unlink( fn ) for fn in os.listdir(cwd) if fn.endswith(".bam.gz") ]
    [ os.unlink( fn ) for fn in os.listdir(cwd) if fn.endswith("sort.gz.bam") ]
    [ os.unlink( fn ) for fn in os.listdir(cwd) if fn.endswith(".bai") ] 

    # remove fastq files
    [ os.unlink( fn ) for fn in os.listdir(cwd) if fn.endswith("fastq.gz") ]
    
    # make wig directory
    if not os.path.exists( cwd + 'wig'):
        os.mkdir( "wig" )
    widDir = cwd + "/wig/"
    [ os.rename( (cwd + fn), (widDir + fn) ) for fn in os.listdir(cwd) if fn.endswith(".wig.gz") ]

    # make HTseq directory
    if not os.path.exists( cwd + 'htseq'):
        os.mkdir( "htseq" )
    htsDir = cwd + "/htseq/"
    [ os.rename( (cwd + fn), (htsDir + fn) ) for fn in os.listdir(cwd) if fn.endswith("_HTseqOutput.txt.gz") ]
    # run HTseqReport.py 
    os.chdir(htsDir)
    cmd = [ '/home/GLBRCORG/mplace/projects/condor/Condor-RNA-Seq-Pipeline/HTseqReport.py' ]
    subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE ).communicate()
    os.chdir(cwd)

    # make log file directory
    if not os.path.exists( cwd + 'log'):
        os.mkdir( "log" )
    logDir = cwd + "/log/"
    [ os.rename( (cwd + fn), (logDir + fn) ) for fn in os.listdir(cwd) if fn.endswith("submit.log") ]
    [ os.rename( (cwd + fn), (logDir + fn) ) for fn in os.listdir(cwd) if fn.endswith("submit.err") ]
    os.rename( (cwd +'pipeline.log'), (logDir + 'pipeline.log') )
    os.rename( (cwd + 'bam2wig.log'), (logDir + 'bam2wig.log') )
    os.rename( (cwd + 'glow.log'), (logDir + 'glow.log') )

    # make a directory for fastqc results
    if not os.path.exists( cwd + 'fastqc'):
        os.mkdir( "fastqc" )
    fastqcDir = cwd + "/fastqc/"
    [ os.rename( (cwd + fn), (fastqcDir + fn) ) for fn in os.listdir(cwd) if fn.endswith("trim_fastqc.zip") ]
    [ os.rename( (cwd + fn), (fastqcDir + fn) ) for fn in os.listdir(cwd) if fn.endswith("trim_fastqc") ]
    
    # delete sam files as they are no longer needed
    [ os.unlink(fn) for fn in os.listdir(cwd) if fn.endswith('.sam.gz')]

def updateGLOW( submitter, rnaDir, wfID, token ):
    """
    Update the submitter's workflow on GLOW with the results of the pipeline.
    Currently only RPKM.results is copied over to GLOW.
    
    Update GLOW with results file location   
    cmd = ['curl', '--cookie',  'cjar', '--data', 'workflow_xml=<workflow workflow_id="59"><datafile><name>RPKM.results</name>\
    <file_path>/mnt/bigdata/processed_data/mplace/RPKM.results</file_path><file_type>Gene-centric Counts</file_type><sub_type>RPKM</sub_type></datafile></workflow>',
           'https://glow-trunk.glbrc.org/upsert_workflow?glow_access_token=c1f3b8d2b8ab8126ad0df366f99a5570d2216b0b' ]
    output = subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()

    Update GLOW with job status (only completion at this point)
    cmd = ['curl', '--cookie', 'cjar', '--data', 'workflow_xml=<workflow workflow_id="130"><completion_status>Complete<
    /completion_status></workflow>',
    'https://glow-trunk.glbrc.org/upsert_workflow?glow_access_token=9a487c01f23738141e7ace6cd1298cbc1c7f720c' ]

    """
    cmd = ['curl', '--cookie',  'cjar', '--data', 'workflow_xml=<workflow workflow_id="' + wfID + '"><datafile><name>RPKM.results</name>\
    <file_path>/mnt/bigdata/processed_data/' + submitter + '/' + rnaDir + '/RPKM.results</file_path><file_type>Gene-centric Counts</file_type><sub_type>RPKM</sub_type></datafile></workflow>',
           'https://glow-trunk.glbrc.org/upsert_workflow?glow_access_token=' + token ]
    output = subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    result1 = output[0].decode('utf-8')
    result2 = output[1].decode('utf-8')

    cmdStatus = ['curl', '--cookie', 'cjar', '--data', 'workflow_xml=<workflow workflow_id="' + wfID +'"><completion_status>Complete</completion_status></workflow>',
    'https://glow-trunk.glbrc.org/upsert_workflow?glow_access_token=' + token ]
    outputStatus = subprocess.Popen( cmdStatus, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    result3 = outputStatus[0].decode('utf-8')
    result4 = outputStatus[1].decode('utf-8')    
    
    with open('glow.log', 'w') as log:
        log.write("GLOW command:")
        log.write(' '.join(cmd))
        log.write("\n")
        log.write(result1)
        log.write("\n")
        log.write(result2)
        log.write("\n")
        log.write(' '.join(cmdStatus))
        log.write("\n")
        log.write(result3)
        log.write("\n")
        log.write(result4)
        log.write("\n")
    log.close()

def bam2wig( bamFile ):
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

def runRPKM( cwd, refer ):
    """
    Run RPKM normalization on all HTSeq files in current working directory.
    
    Calls RPKM.py

    RPKM.py -d <directory> -g REF.gff

    The assumption is that as part of this pipeline HTSeq has been run
    using CDS as the genome feature.
    """
    program = '/home/GLBRCORG/mplace/projects/condor/Condor-RNA-Seq-Pipeline/RPKM.py'
    gff     =  r.ref[refer][2] 
    cmd     = [ program, '-r', refer ]
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
    reference  = sys.argv[1]  # reference genome to use
    submitter  = sys.argv[2]  # users name
    rnaDir     = sys.argv[3]  # rna processing directory name
    workflowID = sys.argv[4]  # workflow ID, used to update GLOW
    token      = sys.argv[5]  # GLOW access token
    currDir    = os.getcwd()
    runRPKM(currDir, reference)
    for file in os.listdir():
        if file.endswith('final.sort.gz.bam'):
            bam2wig(file)
            
    updateGLOW( submitter, rnaDir, workflowID, token )
    cleanUp(currDir, submitter, rnaDir , workflowID, token)
    mail.send("RNA-Seq processing complete")

if __name__ == "__main__":
    main()
 
