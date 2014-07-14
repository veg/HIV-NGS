import os, argparse, re, sys
import subprocess, time, csv
import datetime, shutil, json
#import xlsx
from os.path import join, splitext, isfile

PGMtag = 'TCAGTTCCGATAACGAT'

def run_sff (base_path, results_path):
    print ("Converting sff to FASTQ on %s " % base_path, file = sys.stderr)
    with open (results_path, "w") as out_file:
        try:
            subprocess.check_call (['/usr/local/bin/sff2fastq', '-o', results_path, base_path]) 
        except subprocess.CalledProcessError as e:
            if os.path.exists (results_path):
                return results_path
            print ('ERROR: SFF conversion failed',e,file = sys.stderr)
            return None
        return results_path

def run_qfilt (in_path, in_qual, results_path, status_path):
    print ("Running qfilt on %s saving to %s" % (in_path, results_path), file = sys.stderr)
    with open (results_path, "w") as out_file:
        with open (status_path, "w") as json_file:
            try:
                subprocess.check_call (['/usr/local/bin/qfilt', '-F', in_path , in_qual, '-q', '10', '-l', '50', '-P', '-', '-j'], stdout = out_file, stderr = json_file) 
                
            except subprocess.CalledProcessError as e:
                print ('ERROR: QFILT call failed failed',e,file = sys.stderr)
                return None
            return results_path

def codon_aligner (in_path, out_path, ref = 'data/C_pol.fas'):
    print ("Running bealign on %s " % in_path, file = sys.stderr)
    bam_out = join (out_path, "aligned.bam")
    discards = join(out_path, 'discards.fna')
    try:
        subprocess.check_call (['/usr/bin/bpsh', '6', '/opt/share/python3.3/bealign', '-r', ref, '-e', '0.8', '-m', 'HIV_BETWEEN_F', '-D', discards, '-R', in_path, bam_out ]) 
    except subprocess.CalledProcessError as e:
        print ('ERROR: bealign call failed failed',e,file = sys.stderr)
        return None, None
    return bam_out, discards
    
def bam_to_fasta (in_path, out_path):
    print ("Running bam2msa on %s " % in_path, file = sys.stderr)
    msa_out = join (out_path, "aligned.msa")
    try:
        subprocess.check_call (['/opt/share/python3.3/bam2msa', in_path, msa_out]) 
    except subprocess.CalledProcessError as e:
        print ('ERROR: bam2msa call failed failed',e,file = sys.stderr)
        return None
    return msa_out
    
def collapse_translate_reads (in_path, out_path):
    print ("Collapsing and translating reads %s " % in_path, file = sys.stderr)
    merged_out = join (out_path, "merged.msa")
    merged_out_prot = join (out_path, "merged_prot.msa")
    try:
        subprocess.check_call (['/opt/share/python3.3/seqmerge', in_path, merged_out]) 
        subprocess.check_call (['/opt/share/python3.3/translate', merged_out, merged_out_prot]) 
    except subprocess.CalledProcessError as e:
        print ('ERROR: Collapse/translate call failed failed',e,file = sys.stderr)
        return None
    return (merged_out, merged_out_prot)   
    
def count_collapsed_reads (in_path, out_path):
    print ("Getting protein coverage info for %s " % in_path, file = sys.stderr)
    merged_json= join (out_path, "prot_coverage.json")
    try:
        subprocess.check_call (['/usr/local/bin/seqcoverage', '-o', merged_json, '-t', 'protein', in_path]) 
    except subprocess.CalledProcessError as e:
        print ('ERROR: Protein coverage call failed',e,file = sys.stderr)
        return None
    return merged_json   

def multinomial_filter (in_path, out_path):
    print ("Running multinomial filter on %s " % in_path, file = sys.stderr)
    filtered_out = join (out_path, "filtered.msa")
    json_out = join (out_path, "rates.json")
    try:
        #-t 0.005 -p 0.999999 -f results/Ionxpress020/filtered.msa -j results/Ionxpress020/rates.json results/Ionxpress020/aligned.msa
        subprocess.check_call (['/home/sergei/Projects/MGM/mcmc.jl', '-t', '0.005', '-p', '0.999999', '-f', filtered_out, '-j', json_out, in_path]) 
    except subprocess.CalledProcessError as e:
        print ('ERROR: multinomial filter call failed failed',e,file = sys.stderr)
        return None
    return filtered_out, json_out    

def get_tn93 (in_path, out_path):
    print ("Running tn93 on %s" % in_path, file = sys.stderr)
    try:
        histogram_out = join (out_path, "tn93.json")
        tn93out = join (out_path, "tn93.txt")
        with open (tn93out, "w") as out_file:
            with open (histogram_out, "w") as json_file:
                subprocess.check_call (['/usr/bin/bpsh', '6', '/usr/local/bin/tn93', '-t', str(0.01), '-c', '-q', in_path], stdout = json_file, stderr = out_file) 
        return histogram_out
    except subprocess.CalledProcessError as e:
        print ('ERROR: tn93 call failed',e,file = sys.stderr)
        return None
        

def main (dir, cache_file, results_dir):
                
    if os.path.exists(cache_file):
        with open (cache_file, "r") as fh:
            previous_run_cache = json.load (fh)
    else:
        previous_run_cache = {}
        
    refs = ['HXB2_prrt','HXB2_gag','HXB2_env']
    genes = ['prrt','gag', 'env']
        
    for root, dirs, files in os.walk(dir):
        for file in files:
            name, ext = splitext (file)
            if ext in ('.fna'):
                base_path = join(root, name)
                base_file = join(root, name + ext)

                if base_path not in previous_run_cache:
                    previous_run_cache [base_path] = {'id' : len (previous_run_cache) + 1}
                    
                print ('Working on %s...' % base_file, file = sys.stderr)
                   
                file_results_dir = os.path.join (results_dir, os.path.basename (root))
                if not os.path.exists (file_results_dir):
                    os.makedirs (file_results_dir) 

                if 'in_fasta' not in previous_run_cache [base_path]:
                    previous_run_cache [base_path]['in_fasta'] = base_file
                    previous_run_cache [base_path]['in_qual'] = join(root, name + ".qual")

                if 'in_fasta' in previous_run_cache [base_path] and  previous_run_cache [base_path] ['in_fasta'] is not None:
                    if 'filtered_fastq' not in previous_run_cache [base_path]:
                        previous_run_cache [base_path] ['filtered_fastq'] = run_qfilt (previous_run_cache [base_path] ['in_fasta'], previous_run_cache [base_path] ['in_qual'], join (file_results_dir, 'qfilt.fna'), join (file_results_dir, 'qfilt.json'))

                for i, gene in enumerate(genes): 
                    if gene not in previous_run_cache [base_path]:
                        previous_run_cache [base_path][gene] = {}
                
                    file_results_dir = os.path.join (results_dir, os.path.basename (root), gene)
                    if not os.path.exists (file_results_dir):
                        os.makedirs (file_results_dir) 
                
                    if i == 0:
                        if 'filtered_fastq' in previous_run_cache [base_path] and previous_run_cache [base_path] ['filtered_fastq'] is not None:
                            if 'aligned_bam' not in previous_run_cache [base_path][gene] or previous_run_cache [base_path][gene]['aligned_bam'] is None:
                                previous_run_cache [base_path][gene] ['aligned_bam'], previous_run_cache [base_path][gene] ['discards'] = codon_aligner (previous_run_cache [base_path] ['filtered_fastq'], file_results_dir, refs[i])
                    else:
                        prev_gene = genes[i-1]
                        if 'discards' in previous_run_cache [base_path][prev_gene] and previous_run_cache [base_path][prev_gene]['discards'] is not None:
                            if 'aligned_bam' not in previous_run_cache [base_path][gene] or previous_run_cache [base_path][gene]['aligned_bam'] is None:
                                previous_run_cache [base_path][gene] ['aligned_bam'], previous_run_cache [base_path][gene] ['discards'] = codon_aligner (previous_run_cache [base_path][prev_gene]['discards'], file_results_dir, refs[i])

                    if 'aligned_bam' in previous_run_cache [base_path][gene] and previous_run_cache [base_path][gene] ['aligned_bam'] is not None:
                        if 'aligned_msa' not in previous_run_cache [base_path][gene]:
                            previous_run_cache [base_path][gene] ['aligned_msa'] = bam_to_fasta (previous_run_cache [base_path][gene] ['aligned_bam'], file_results_dir)

                    if 'aligned_msa' in previous_run_cache [base_path][gene] and previous_run_cache [base_path][gene] ['aligned_msa'] is not None:
                        if 'filtered_msa' not in previous_run_cache [base_path][gene] or 'json_rates' not in previous_run_cache [base_path][gene]:
                            previous_run_cache [base_path][gene] ['filtered_msa'], previous_run_cache [base_path][gene] ['json_rates']  = multinomial_filter (previous_run_cache [base_path][gene] ['aligned_msa'], file_results_dir)
        
                    if 'filtered_msa' in previous_run_cache [base_path][gene] and previous_run_cache [base_path][gene] ['filtered_msa'] is not None:
                        if 'merged_msa' not in previous_run_cache [base_path][gene] or 'merged_msa_prot' not in previous_run_cache [base_path][gene]:
                            previous_run_cache [base_path][gene] ['merged_msa'], previous_run_cache [base_path][gene] ['merged_msa_prot']  = collapse_translate_reads (previous_run_cache [base_path][gene] ['filtered_msa'], file_results_dir)
 
                    if 'merged_msa_prot' in previous_run_cache [base_path][gene] and previous_run_cache [base_path][gene] ['merged_msa_prot'] is not None:
                        if 'merged_json' not in previous_run_cache [base_path][gene]:
                            previous_run_cache [base_path][gene] ['merged_json'] = count_collapsed_reads (previous_run_cache [base_path][gene] ['merged_msa_prot'], file_results_dir)

                    if 'merged_msa' in previous_run_cache [base_path][gene] and previous_run_cache [base_path][gene] ['merged_msa'] is not None:
                        if 'tn93_json' not in previous_run_cache [base_path][gene]:
                            previous_run_cache [base_path][gene] ['tn93_json'] = get_tn93 (previous_run_cache [base_path][gene] ['merged_msa'], file_results_dir)

                dthandler = lambda obj: obj.isoformat() if isinstance(obj, datetime.datetime) else None  
                with open (cache_file, "w") as fh:
                    json.dump (previous_run_cache, fh,default=dthandler, sort_keys=True, indent=4)

    return 0

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='scan the directory of NGS files and process them'
    )
    parser.add_argument(
        '-i', '--input',
        metavar='DIR',
        type=str,
        help='the directory to scan',
        required = True,
    )
    parser.add_argument(
        '-c', '--cache',
        metavar='JSON',
        type=str,
        help='the file which contains the .json cache file',
        required = True,
    )
    parser.add_argument(
        '-r', '--results',
        metavar='RESULTS',
        type=str,
        help='the directory where analysis results will be written',
        required = True,
    )    
    
    args = None
    retcode = -1
    args = parser.parse_args()
    
    if not os.path.exists (args.results):
        os.mkdir (args.results) 
    
    retcode = main(args.input, args.cache, args.results)
    
    sys.exit(retcode)


