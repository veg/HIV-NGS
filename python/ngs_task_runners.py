import os, argparse, re, sys
import subprocess, time, csv
import datetime, shutil, json, operator
#import xlsx
from os.path import join, splitext, isfile
from Bio import SeqIO
import threading
import queue
import ngs_task_runners


check_file_paths = True
check_file_paths_diversity = True

def set_cache_flags (f1, f2):
    global check_file_paths
    global check_file_paths_diversity
    check_file_paths = f1
    check_file_paths_diversity = f2
    
    
def codon_aligner (in_path, out_path, node, ref = 'HXB2_prrt'):
    print ("Running bealign on %s (node %d)" % (in_path, node), file = sys.stderr)
    bam_out = join (out_path, "aligned.bam")
    discards = join(out_path, 'discards.fna')

    if check_file_paths and os.path.exists (bam_out) and os.path.getsize ( bam_out ) and  os.path.exists (discards):
        return bam_out, discards
          
    command = ['/usr/bin/bpsh', str(node), '/opt/share/pytools/utils/python3.3/bealign', '-r', ref, '-e', '0.5', '-m', 'HIV_BETWEEN_F', '-D', discards, '-R', in_path, bam_out ]      
    try:
        subprocess.check_call (command, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL) 
    except subprocess.CalledProcessError as e:
        print ('ERROR: bealign call failed failed',' '.join (command),file = sys.stderr)
        return None, None
    return bam_out, discards if isfile (discards) else None 
    
def bam_to_fasta (in_path, out_path):
    print ("Running bam2msa on %s " % in_path, file = sys.stderr)
    msa_out = join (out_path, "aligned.msa")
    sorted_bam = join (out_path, "kludge")
    if check_file_paths and os.path.exists (msa_out) and os.path.getsize ( msa_out ) > 0:
        return msa_out
    try:
        print ("Running samtools sort on %s to %s" % (in_path,sorted_bam), file = sys.stderr)
        subprocess.check_call (['/usr/local/bin/samtools', 'sort', in_path, sorted_bam],stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL) 
        shutil.move (sorted_bam + ".bam", in_path)
        subprocess.check_call (['/opt/share/python3.3/bam2msa', in_path, msa_out],stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL) 
    except subprocess.CalledProcessError as e:
        print ('ERROR: bam2msa call failed failed',e,file = sys.stderr)
        return None
    return msa_out
    

def get_tn93 (in_path, out_path, node ):
    print ("Running tn93 on %s (node %d)" % (in_path, node), file = sys.stderr)
    try:
        histogram_out = join (out_path, "tn93.json")
        tn93out = join (out_path, "tn93.txt")

        if check_file_paths and os.path.exists (histogram_out) and os.path.exists (tn93out):
            return histogram_out

        with open (tn93out, "w") as json_file:
            with open (histogram_out, "w") as out_file:
                subprocess.check_call (['/usr/bin/bpsh', str(node), '/usr/local/bin/tn93', '-t', str(0.01), '-c', '-q', in_path], stdout = json_file, stderr = out_file) 
        return histogram_out
    except subprocess.CalledProcessError as e:
        print ('ERROR: tn93 call failed',e,file = sys.stderr)
        return None




