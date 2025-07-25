#!/usr/bin/env python

# 2018 - Centre de Regulacio Genomica (CRG) - All Rights Reserved

# Writen by Miravet-Verde, Samuel
# Last updated = 12/10/2018 (MM/DD/YYYY)

# Command line example:
# python /software/ls/fastqins/fastqins_web.py -i1 /users/lserrano/smiravet/ins_results/web_tests/test_read1.fastq.gz -i2 /users/lserrano/smiravet/ins_results/web_tests/test_read2.fastq.gz -t TACGGACTTTATC -g /users/lserrano/www/reference_genome/mycoplasma_pneumoniae_m129.gbk -o /users/lserrano/smiravet/ins_results/web_tests/your_test/ -p 1 -v -r 0

import glob
import json
import os.path
import sys, os
import argparse
import subprocess
from Bio import SeqIO

# --------------------------------------------------

# Set some program accesses
# Update to point to the path with your fastuniq, bowtie and samtools installation
softdir  = '/software/ls/'
sys.path.insert(0, softdir+"fastqins")
from socketIO_client import SocketIO, LoggingNamespace

fastuniq = softdir+'FastUniq/source/fastuniq'
bowtie2  = softdir+'bowtie2-2.2.9/bowtie2'
samtools = softdir+'/samtools-1.4/samtools'

####
# FUNCTIONS
####

def check_everything(lista):
    for fname in lista:
        if not os.path.isfile(fname):
            print("ERROR: no input file exists with path(s)", '\n'.join([str(i) for i in lista]))
            raise SystemExit

def revcomp(seq):
    comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    return ''.join([comp[i] for i in seq])[::-1]

def testfasta_and_length(genome, outdir):
    genome_id = genome.split('/')[-1].split('.')[0]
    genome_ex = genome.split('/')[-1].split('.')[-1]

    # Determine the file type:
    if genome_ex in ['gb', 'gbk', 'genbank']:
        tipo = 'genbank'
    elif genome_ex in ['fa', 'fna', 'fasta', 'fast']:
        tipo = 'fasta'
    else:
        print("GENOME FORMAT NOT SUPPORTED")
        raise SystemExit

    # Write and return
    handle = open(genome, 'rU')
    for record in SeqIO.parse(handle, tipo):
        genome_size = str(len(record.seq))
        if tipo=='fasta':
            return genome, genome_size
        elif tipo=='genbank':
            fout = outdir+genome_id+'.fna'
            handleout = open(fout, 'w')
            handleout.write(">{}\n{}".format(record.id, record.seq))
            handleout.close()
            return fout, genome_size
    handle.close()

def project_identifier(fil, outputfolder):
    """
    Required to discriminate between @
    in the identifier and in the quality string
    """
    out = outputfolder+'id.txt'
    cmd='head -1 '+fil+' > '+out
    os.system(cmd)
    with open(out, 'rU') as fi:
        for line in fi:
            return line.split(':')[0]

def extract_insertions(out, read_length, genome_size):
    """ core function to call insertions """

    # Filter by read length: only take reads that are shorter than usual: they have the transposon trimmed!
    # The position is equal to the first base mapped before the transposon starts 
    # 'XS:' --> remove reads with multiple mapping

    # Forward strand
    cmd = "awk 'length($10) <"+read_length+" || $1 ~ /^@/' "+out+"_forward_filt_paired.sam | grep -v '^@' | grep -v 'XS:' | awk '{print $4}' > "+out+"_forward.ins"
    os.system(cmd)
    # Reverse strand
    cmd = "awk 'length($10) <"+read_length+" || $1 ~ /^@/' "+out+"_reverse_filt_paired.sam | grep -v '^@' | grep -v 'XS:' | awk '{if($4+length($10)-1<="+genome_size+") {print $4+length($10)-1} else if ($4+length($10)-1>"+genome_size+") {print $4+length($10)-1-"+genome_size+"}}' > "+out+"_reverse.ins"
    os.system(cmd)

    # Merge and extract info:
    cmd = "cat "+out+"_reverse.ins "+out+"_forward.ins | sort | uniq -c | awk '{print $2, $1}' > "+out+".qins"
    os.system(cmd)

    #@@ NEW # separate by strand
    cmd = "cat "+out+"_reverse.ins | sort | uniq -c | awk '{print $2, $1}' > "+out+"_reverse.qins"
    os.system(cmd)
    cmd = "cat "+out+"_forward.ins | sort | uniq -c | awk '{print $2, $1}' > "+out+"_forward.qins"
    os.system(cmd)


def extract_identifier_position(out, read_length, genome_size):
    # Extracting barcode information
    cmd = "awk 'length($10) <"+read_length+" || $1 ~ /^@/' "+out+"_forward_filt_paired.sam | grep -v '^@' | grep -v 'XS:' | awk '{print $1, $4}' > "+out+"_forward.ides"
    os.system(cmd)
    cmd = "awk 'length($10) <"+read_length+" || $1 ~ /^@/' "+out+"_reverse_filt_paired.sam | grep -v '^@' | grep -v 'XS:' | awk '{if($4+length($10)-1<="+genome_size+") {print $1, $4+length($10)-1} else if ($4+length($10)-1>"+genome_size+") {print $1, $4+length($10)-1-"+genome_size+"}}' > "+out+"_reverse.ides"
    os.system(cmd)
    cmd = "cat "+out+"_reverse.ides "+out+"_forward.ides > "+out+".ides"
    os.system(cmd)

def finish_log_file(fastq_dir, fastq_filt, tn_seq, out, genome_size):

    # Compute stats:
    cmd = "awk '{s++}END{print s/4}' "+fastq_dir
    nr_reads = float(subprocess.check_output(cmd, shell=True).strip())
    if fastq_filt:
        cmd = "awk '{s++}END{print s/4}' "+fastq_filt
        nr_remain = float(subprocess.check_output(cmd, shell=True).strip())
        filt_line = str(int(nr_remain))+" ("+str(round(100.0*nr_remain/nr_reads, 2))+"%) passed the PCR duplicates filter; of these:"

        cmd = "grep "+tn_seq+" "+fastq_filt+" | wc -l"
        nr_IR = float(subprocess.check_output(cmd, shell=True).strip())
    else:
        nr_remain = nr_reads
        filt_line = "PCR duplicate removal not selected. "+str(int(nr_reads))+" considered; of these:"

        cmd = "grep "+tn_seq+" "+fastq_dir+" | wc -l"
        nr_IR = float(subprocess.check_output(cmd, shell=True).strip())

    cmd = "cat "+out+"_reverse.ins "+out+"_forward.ins | wc -l"
    nr_Tn = float(subprocess.check_output(cmd, shell=True).strip())
    cmd = "cat "+out+".qins | wc -l"
    nr_pos = int(subprocess.check_output(cmd, shell=True).strip())
    if nr_pos!=0:
        cmd = "awk '{sum+=$2} END{print sum}' "+out+".qins"
        nr_post_reads = int(subprocess.check_output(cmd, shell=True).strip())
    else:
        nr_post_reads = 0
    # Write log
    separator = '\n-----------------\n\n'
    with open(out+'.log', 'a') as logfile:
        logfile.write(separator+'FASTQINS INFORMATION:\n\n')
        text = "{0} reads provided; of these:\n  ".format(int(nr_reads))+filt_line
        text+= "\n    {0} ({1}%) presented the IR:{2} sequence; of these:\n      ".format(int(nr_IR), round(100.0*nr_IR/nr_remain, 2), tn_seq)
        text+= "{0} having IR ({1}%) were mapped unambiguously; in total:\n        ".format(int(nr_Tn), round(100.0*nr_Tn/nr_IR, 2))
        text+= "{0} insertion positions can be extracted from {1} reads (this number should be the same than previous line)\n".format(nr_pos, nr_post_reads)
        text+= "{0}% of the fastq was informative.".format(round(100.0*nr_Tn/nr_reads, 2))
        text+= "{0}% coverage (number of insertions per genome base).".format(round(100*int(nr_pos)/int(genome_size)))
        logfile.write(text)


def fastqins(read1        , read2   , Tn_seq      , genome    ,
             output_folder, pcr_dup , barcode_size, mismatches,
             extension    , rm_inter, verbose    ):

    check_everything([read1, read2])

    Tn_seq = Tn_seq.upper()

    if verbose:
        print('read1:', read1)
        print('read2:', read2)
        print('TnSeq:', Tn_seq)
        rcTn_seq = revcomp(Tn_seq)
        print('RC-TnSeq:', rcTn_seq)
        print('genome:', genome)
        print('outFolder:', output_folder)

    # Create the directory and the environment
    if output_folder[-1]!='/':
        output_folder+='/'
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    output_folder += extension+read1.split('/')[-1].split('.')[0]+'_intermediate_files/'
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # If genome in genbank format, create a fasta file
    genome, genome_size = testfasta_and_length(genome, output_folder)

    # Create log file:
    pipeline_name = 'pipeline_fastqins'
    pipelineDoc = pipeline_name + """
    Pipeline to analyze Tn-seq data.
    Python >=2.7 script. |
    Author: Miravet-Verde, Samuel |
    Last updated: 2018.12.10 |
    Affiliation: Center for Genomic Regulation, Luis Serrano's lab |
    email: samuelmiver@gmail.com; samuel.miravet@crg.eu |\n
    """
    separator = '\n-----------------\n\n'
    with open(output_folder+read1.split('/')[-1].split('.')[0]+'.log', 'w') as logfile:
        logfile.write(pipelineDoc)
        logfile.write('Execution info\n-----------------\nversion:FASTQINS v1.0'+separator+'GENERAL INFORMATION:\n\nR1:{0}\nR2:{1}\nTnSeq:{2}\nGenome:{3}\nGenome length:{4}\nOutDir:{5}\nPCR dup removal:{6}\nMismatches:{7}'.format(read1, read2, Tn_seq, genome, genome_size, output_folder, str(pcr_dup), str(mismatches))+separator+'BOWTIE2 ALIGN. INFORMATION:\n\n')

    midfile1 = output_folder+extension+read1.split('/')[-1]
    midfile2 = output_folder+extension+read2.split('/')[-1]

    # Decompress the file if compressed
    marker = False
    if read1.endswith('gz'):
        cmd = 'gunzip -c '+read1+' > '+midfile1.replace('.gz', '')
        if verbose:
            print('decompressing read1 ... ')
        _ = os.system(cmd)
        midfile1 = midfile1.replace('.gz', '')
        marker = True

    if read2.endswith('gz'):
        cmd = 'gunzip -c '+read2+' > '+midfile2.replace('.gz', '')
        if verbose:
            print('decompressing read2 ... ')
        _ = os.system(cmd)
        midfile2 = midfile2.replace('.gz', '')
        read2_backup = str(midfile2)
    else:
        read2_backup = str(read2)

    # Get the read length
    cmd = 'head -2 '+read2_backup+' | tail -1'
    read_length = str(len(subprocess.check_output(cmd, shell=True))-1)

    # Remove duplicates
    read2_filt = False
    if pcr_dup:
        print('removing duplicates')
        with open(output_folder+'lis.ls', 'w') as fo:
            if marker:
                fo.write(midfile1+'\n'+midfile2)
            else:
                fo.write(read1+'\n'+read2)
        if verbose:
            print('Removing PCR duplicates...')
        cmd = fastuniq+' -i '+output_folder+'lis.ls -t q -o '+midfile1+'.filt -p '+midfile2+'.filt'
        _ = os.system(cmd)
        print "Command output: " + output

        if verbose:
            print('PCR Duplicates removed.')
        midfile1 += '.filt'
        midfile2 += '.filt'
        read2_filt = str(midfile2)


    #@@NEW
    # IDENTIFIERS WITH TN

    # This block is old we added the filtering grep at the beginning
    if marker or pcr_dup:
        cmd = "grep -B1 '"+Tn_seq+"' "+midfile2+" | grep '@' > "+midfile2+'.ides'
        os.system(cmd)

        cmd0 = "java -jar /users/lserrano/smiravet/soft/Trimmomatic-0.38/trimmomatic-0.38.jar SE -phred33 "+midfile2+" "+midfile2+".trim ILLUMINACLIP:/users/lserrano/smiravet/soft/Trimmomatic-0.38/adapters/TruSeq3-SE.fa:2:30:10"
        midfile2 += '.trim'
    else:
        cmd = "grep -B1 '"+Tn_seq+"' "+read2+" | grep '@' > "+midfile2+'.ides'
        os.system(cmd)

        cmd0 = "java -jar /users/lserrano/smiravet/soft/Trimmomatic-0.38/trimmomatic-0.38.jar SE -phred33 "+read2+" "+midfile2+".trim ILLUMINACLIP:/users/lserrano/smiravet/soft/Trimmomatic-0.38/adapters/TruSeq3-SE.fa:2:30:10"
        midfile2 += '.trim'
    if verbose:
        print('Calling the transposon and trimming it ... ')

    _ = os.system(cmd0)
    cmd = "grep -F -w -A3 -f "+midfile2.replace('.trim', '.ides')+" "+midfile2+" | grep -v '^--$' | awk 'NR%4==2{s=match($0,/"+Tn_seq+"/)+RLENGTH} NR%4~/[02]/{$0=substr($0,s)} 1' > "+midfile2+".tn"
    os.system(cmd)
    midfile2 += '.tn'

    # Barcode processement if present
    if barcode_size:
        if verbose:
            print('Calling the barcode with size '+str(barcode_size)+' and trimming it ...')
        # Generate a fasta with the identifiers and the associated barcode
        # Trim the barcode in the read1 file (when line does not start with @IDE or +
        if marker or pcr_dup:
            ide = project_identifier(midfile1, output_folder)
            cmd1 = "sed -n '1~4s/^"+str(ide)+"/>/p;2~4p' "+midfile1+" | awk '{if($0~/^>/) {print $0} else {print substr($0,0,"+str(barcode_size)+")}}' > "+midfile1+"_barcodes.fa"
            cmd2 = "awk '{if($0 ~/^"+str(ide)+"/ || $0 ~/+/) {print $0} else {print substr($0,"+str(int(barcode_size)+1)+")}}' "+midfile1+" > "+midfile1+".bc"
        else:
            ide = project_identifier(read1, output_folder)
            cmd1 = "sed -n '1~4s/^"+str(ide)+"/>/p;2~4p' "+read1+" | awk '{if($0~/^>/) {print $0} else {print substr($0,0,"+str(barcode_size)+")}}' > "+midfile1+"_barcodes.fa"
            cmd2 = "awk '{if($0 ~/^"+str(ide)+"/ || $0 ~/+/) {print $0} else {print substr($0,"+str(int(barcode_size)+1)+")}}' "+read1+" > "+midfile1+".bc"
        os.system(cmd1)
        os.system(cmd2)
        midfile1 += '.bc'

    # Index the genome
    if verbose:
        print('Indexing the genome '+genome+' ...')
    cmd = bowtie2+'-build -f '+genome+' '+output_folder+'genome_bowtie'
    os.system(cmd)

    # Map the reads
    out = output_folder+midfile1.split('/')[-1].split('.')[0]
    # -p : threads
    # -x : genome
    # -1 : read1
    # -2 : read2
    # -S : sam output
    cmd = bowtie2+' -x '+output_folder+'genome_bowtie -N '+str(mismatches)+' -U '+midfile2+' -S '+out+'.sam 2>>'+out+'.log'
    if verbose:
        print('Mapping ...')
    os.system(cmd)

    # Transform to bam
    # -b : bam output
    # -S : sam input
    # -o : output name
    cmd = samtools+' view -b -S -o '+out+'.bam '+out+'.sam'
    if verbose:
        print('Creating bam ...')
    os.system(cmd)

    # Filter paired end
    # -F 0x04 : filtering unmapped reads
    # -f 0x02 : required to map PE
    # -q 30 : minimum alignment quality
    # cmd = samtools+' view -F 0x04 -f 0x02 -q 30 -b '+out+'.bam > '+out+'.paired_mapped'
    cmd = samtools+' view -q 10 -F 0x10 '+out+'.bam > '+out+'_forward_filt_paired.sam'   # Analog to -F 16, 0x10 is 16 in decimal value
    os.system(cmd)
    cmd = samtools+' view -q 10 -f 0x10 '+out+'.bam > '+out+'_reverse_filt_paired.sam'
    os.system(cmd)
    if verbose:
        print('Filtering unmmaped and separating forward/reverse...')

    # Take unpaired but mapped reads
    # -F 0x04 : filtering unmapped reads
    # -F 0x02 : filtering PE reads
    # cmd = samtools+' view -F 0x04 -F 0x02 -b '+out+'.bam > '+out+'.unpaired_mapped'
    # print 'Extracting mapped that are non-paired...'
    # os.system(cmd)
    # Back to sam the 2 files
    # cmd = samtools+' view -h -o '+out+'_filt_paired.sam '+out+'.paired_mapped'
    # print 'Back to sam the paired mapped reads ...'
    # os.system(cmd)
    # cmd = samtools+' view -h -o '+out+'_filt_single.sam '+out+'.unpaired_mapped'
    # print 'Back to sam the singletons mapped ...'
    # os.system(cmd)

    # Extract insertions and relation with identifiers
    if verbose:
        print('Selecting reads with the transposon from paired reads N<'+str(read_length)+'...')
    extract_insertions(out, read_length, genome_size)
    if barcode_size:
        if verbose:
            print('Extracting barcode information...')
        extract_identifier_position(out, read_length, genome_size)
        cmd = 'cp '+output_folder+'*.ides '+output_folder.replace(extension+read1.split('/')[-1].split('.')[0]+'_intermediate_files/', '')
        os.system(cmd)

    # Finish the log file
    finish_log_file(read2_backup, read2_filt, Tn_seq, out, genome_size)
    # Put important files in its place
    cmd = 'cp '+output_folder+'*.qins '+output_folder.replace(extension+read1.split('/')[-1].split('.')[0]+'_intermediate_files/', '')
    os.system(cmd)
    cmd = 'cp '+output_folder+'*.log '+output_folder.replace(extension+read1.split('/')[-1].split('.')[0]+'_intermediate_files/', '')
    os.system(cmd)

    if rm_inter and '_intermediate_files' in output_folder:
         cmd = 'rm '+output_folder+'*'
         os.system(cmd)
         cmd = 'rmdir '+output_folder
         os.system(cmd)

         if verbose:
             print('Intermediate files removed.\n')

    # Say goodbye
    if verbose:
        print('Insertions extracted to '+out+".qins file\n\nEnjoy your day and remember:\n\nSTAND AND BE TRUE :D!\n")


# FOR WEB SERVER
def communicate_with_socket(analysisId, status, port=50001):
    with SocketIO('dbspipes.crg.es', port, LoggingNamespace) as socketIO:
        # Send message using web socket to the web server DBSpipes
        print("Sending message to web server via websocket, analysisId="+str(analysisId)+ " and status="+str(status))
        data = {"internal_id": analysisId, "status": status}
        socketIO.emit('on_update', json.dumps(data))
        # Listen
        socketIO.wait(seconds=1)

def write_jobid_files(OUTFILES, analysisId, bc):
    # Check that final output files exist
    # task_path is the last task's output directory
    if bc!=0:
        n = 3
        indexes, extensions = [0,1,2], ['.qins', '.log', '.ides']
    else:
        n = 2
        indexes, extensions = [0,1], ['.qins', '.log']

    if len(OUTFILES)!=n:
        status = -1
    else:
        for index, extension in zip(indexes, extensions):
            # We just check that the coverage qins file exists and is not empty
            if OUTFILES[index].endswith(extension) and os.path.getsize(OUTFILES[index]) > 0:
                status = 2
            else:
                status = -1
    communicate_with_socket(analysisId, status)

# --------------------------------------------------------------

#####
# PARSER
#####

parser = argparse.ArgumentParser(description = "FastQins extracts insertion positions from a fastq file. Please cite ''")

parser.add_argument('-i1', '--read1',
                    dest="read1",
                    required=True,
                    type=str,
                    help="Read 1, no transposon in this read.")

parser.add_argument('-i2', '--read2',
                    dest="read2",
                    required=True,
                    type=str,
                    help="Read 2, transposon IR in this read.")

parser.add_argument('-t', '--Tn_seq',
                    dest="Tn_seq",
                    default="TACGGACTTTATC",
                    type=str,
                    help="Inverted Repeat Transposon sequence to trim. Default= TACGGACTTTATC")

parser.add_argument('-g', '--genome',
                    dest="genome",
                    default="/users/lserrano/www/reference_genome/mycoplasma_pneumoniae_m129.gbk",
                    type=str,
                    help="Genome of reference to map insertions. Default= /users/lserrano/www/reference_genome/mycoplasma_pneumoniae_m129.gbk")

parser.add_argument('-o', '--output_folder',
                    dest="output_folder",
                    default="/users/lserrano/smiravet/ins_results/web_tests/",
                    type=str,
                    help="Output directory to write the files.")

parser.add_argument('-p', '--pcr_dup',
                    dest="pcr_dup",
                    default=1,
                    type=int,
                    help="Either to remove or not the pcr dup. Default=True. Not recommended in highly selected passages.")

parser.add_argument('-b', '--barcode_size',
                    dest="barcode_size",
                    default=0,
                    type=int,
                    help="Extract barcodes or not. [NOT TESTED IN EVERY BARCODE CONDITION]")

parser.add_argument('-m', '--mismatches',
                    dest="mismatches",
                    default=0,
                    type=int,
                    help="Accepted mismacthes during alignment. Default=0.")

parser.add_argument('-e', '--extension',
                    dest="extension",
                    default='',
                    type=str,
                    help="Extension to add to the name of the file. Default='' (nothing).")

parser.add_argument('-r', '--rm_inter',
                    dest="rm_inter",
                    default=1,
                    type=int,
                    help="Remove intermediate files. By default intermediate files are removed.")

parser.add_argument("-v", "--verbose",
                    dest="verbose",
                    action="store_true",
                    help="increase output verbosity")

# TO RUN IN WEB
parser.add_argument('-id', '--analysisId',
                    dest="analysisId",
                    action="store",
                    type=str,
                    help="Custom analysis ID for the cluster.")

parser.add_argument('-w', '--run_on_cluster',
                    dest="run_on_cluster",
                    action="store_true",
                    help="Flag to run on cluster.")

# LOAD:
options = parser.parse_args()

# --------------------------------------------------------------
####
# MAIN RUN
####


# If in cluster
if options.run_on_cluster:
    check_everything([options.read1, options.read2])
    # Check directory structure
    if options.output_folder[-1]!='/':
        new_output = options.output_folder+'/'
    if not os.path.exists(new_output):
        os.makedirs(new_output)

    # Write bash script
    bash_script = '#!/bin/sh\n\n'
    bash_script += 'python /software/ls/fastqins/fastqins_web.py -i1 '+str(options.read1)+' -i2 '+str(options.read2)
    bash_script += ' -t '+str(options.Tn_seq)+' -g '+str(options.genome)+' -o '+str(new_output)+' -p '+str(options.pcr_dup)+' -v'
    bash_script += ' -id '+str(options.analysisId)+' -r '+str(options.rm_inter)
    fname = new_output+"qsub_fastqins_"+str(options.analysisId)+".sh"
    with open(fname, "w") as fo:
        fo.write(bash_script)

    # Submit to queue
    if ' -w ' not in bash_script:
        cmd = 'qsub -q long-sl7 -N qsub_fastqins_'+str(options.analysisId)+' -o '+new_output+' -e '+new_output+' '+fname
        # cmd = 'qsub -q long-sl7 -l virtual_free=24G,h_rt=12:00:00 -pe smp 8 -N qsub_fastqins_'+str(options.analysisId)+' -o '+new_output+' -e '+new_output+' '+fname
        # cmd = 'qsub -q long-sl7 -l virtual_free=48G,h_rt=48:00:00 -pe smp 8 -N qsub_fastqins_'+str(options.analysisId)+' -o '+new_output+' -e '+new_output+' '+fname
        os.system(cmd)
    else:
        print("ERROR: double web call")
        raise SystemExit
else:
    # This try-except is required to communicate with socket if killed by keyboard.
    try:
        # The job entered the queue at this step, send status
        if options.analysisId:
            status = 1
            communicate_with_socket(options.analysisId, status)

        # Run
        fastqins(read1=options.read1  , read2=options.read2,
                 Tn_seq=options.Tn_seq,
                 genome=options.genome,
                 output_folder=options.output_folder,
                 pcr_dup=options.pcr_dup, barcode_size=options.barcode_size,
                 mismatches=options.mismatches, extension=options.extension,
                 rm_inter=options.rm_inter, verbose=options.verbose)
        # Check output and send messages if in cluster
        if options.analysisId:
            handle = options.output_folder
            if handle[-1]!='/':
                handle+='/'
            OUTFILES =  glob.glob(handle+'*.qins')
            OUTFILES += glob.glob(OUTFILES[0].replace('.qins','.log'))
            if options.barcode_size:
                OUTFILES += glob.glob(handle+'*.ides')
            print OUTFILES
            write_jobid_files(OUTFILES, analysisId=options.analysisId, bc=options.barcode_size)
    except:
        if options.analysisId:
            communicate_with_socket(options.analysisId, -1)

# 2018 - Centre de Regulacio Genomica (CRG) - All Rights Reserved
