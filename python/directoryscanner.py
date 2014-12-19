#------------------------------------------------------------------------------#
#                             directoryscanner.py                              #
#------------------------------------------------------------------------------#




# First, import dependencies

import os, argparse, re, sys
import subprocess, time, csv
import datetime, shutil, json, operator, hashlib
#import xlsx
from os.path import join, splitext, isfile
from Bio import SeqIO
import threading
import queue
import ngs_task_runners as ntr
from itertools import product, combinations




# Declare global variables

global previous_run_cache
global threading_lock
global task_queue
global force_diversity_estimation



#------------------------------------------------------------------------------#
#                          Initialize some variables                           #
#------------------------------------------------------------------------------#

# check_file_paths appears later in the multinomial_filter function.
# check_file_paths_diversity appears later in the functions 
# collapse_translate_reads, count_collapsed_reads, extract_diagnostic_region, 
# and collapse_diagnostic_region. path_to_this_file is the path to the
# directory containing directoryscanner.py.

check_file_paths = True
check_file_paths_diversity = True
path_to_this_file = os.path.dirname(os.path.realpath(__file__))

# applying ngs_task_runners.set_cache_flags sets check_file_paths and 
# check_file_paths_diversity as global variables.

ntr.set_cache_flags (check_file_paths, check_file_paths_diversity)

# known_refs is a list of paths to reference files for the genes indicated in 
# the known_genes list, below. These references will be used for aligning 
# codons later. Reference sequences are accessed during handle_a_gene, for 
# alignment and for updating the analysis cache. known_refs and known_genes are
# used early in the main script to check that all of the input genes are 
# included in known_genes. genes is an empty list to be filled by the argument
# parser.

known_refs  = [os.path.normpath(os.path.join (path_to_this_file, k)) for k in ['../data/rt.fas','../data/gag_p24.fas','../data/env_C2V5.fas', '../data/pr.fas', '../data/hcv1bcore.fas']]
known_genes = ['rt','gag', 'env', 'pr', 'hcv1bcore']
genes       = []

# cache_file is a file storing cached information from update_global_record, 
# compartmentalization_handler, and the main script. It is initialized here as 
# an empty string to be filled by the argument parser. dthandler is used during
# json dumps to the cache file. nodes_to_run_on is a list of nodes to run on,
# initialized here as an empty list to be filled by the argument parser.

cache_file = ""
dthandler = lambda obj: obj.isoformat() if isinstance(obj, datetime.datetime) else None  
nodes_to_run_on = []

# spans is initialized here as an empty list to be filled by the argument 
# parser. handle_a_gene iterates on spans when calling
# extract_diagnostic_region. window is used as a constant by
# extract_and_collapse_well_covered_region as an argument to 
# extract_diagnostic_region. stride is used as a constant when updating spans.

spans = []
window   = 210
stride   = 30

# force_diversity_estimation is a toggle for various steps in handling the 
# results list. force_qfilt_rerun causes qfilt to be run automatically, 
# regardless of whether filtered_fastq is in NGS_run_cache.

force_diversity_estimation = False
force_qfilt_rerun           = False




#------------------------------------------------------------------------------#
#                              Define Functions                                #
#------------------------------------------------------------------------------#

#------------------------------- Task Runners ---------------------------------#

### run_sff ###

# Convert sff to FASTQ. This function takes a path to an input .sff file and a 
# path to an output .fastq file as arguments. At present it is not used in the
# main script.

def run_sff (base_path, results_path):
    print ("Converting sff to FASTQ on %s " % base_path, file = sys.stderr)
    
    # Open and write to the output .fastq file.
    
    with open (results_path, "w") as out_file:
        
        # Call the sff2fastq utility on the input .sff file and write to the
        #output .fasta.
        
        try:
            subprocess.check_call (['/usr/local/bin/sff2fastq', '-o', results_path, base_path]) 
            
        # If an error occurs, say so and return the path to the output file if 
        # possible.
        
        except subprocess.CalledProcessError as e:
            if os.path.exists (results_path):
                return results_path
            print ('ERROR: SFF conversion failed',e,file = sys.stderr)
            return None
            
        # Return the path to the output file.
            
        return results_path



        
### run_qfilt ###

# Apply qfilt to remove or split reads with low quality scores. This function
# takes as its arguments either paths to an input .fna and .qual pair or a 
# path to a single input .fastq file, along with a path to an output file
# and an error file.

def run_qfilt (in_path, in_qual, results_path, status_path):
    print ("Running qfilt on %s saving to %s" % (in_path, results_path), file = sys.stderr)
    
    # Open and write to the output and error files.
    
    with open (results_path, "w") as out_file:
        with open (status_path, "w") as json_file:
            
            # If an input .fna and .qual pair is given, apply the qfilt utility with
            # the -F flag. If an input .fastq is given, apply qfilt with the -Q flag.
            
            try:
                if in_qual is not None:
                    subprocess.check_call (['/usr/local/bin/qfilt', '-F', in_path , in_qual, '-q', '15', '-l', '50', '-P', '-', '-R', '8', '-j'], stdout = out_file, stderr = json_file) 
                else:
                    subprocess.check_call (['/usr/local/bin/qfilt', '-Q', in_path , '-q', '15', '-l', '50', '-P', '-', '-j', '-R', '8'], stdout = out_file, stderr = json_file) 
            
            # If an error occurs, say so and return the path to the output file if
            # possible.
            
            except subprocess.CalledProcessError as e:
                print ('ERROR: QFILT call failed failed',e,file = sys.stderr)
                try:
                    if os.path.getsize (results_path) > 2**16:
                        return results_path
                except:
                    pass
                    
                return None
            
            # Return the path to the output file.    
            
            return results_path




### collapse_translate_reads ###

# Collapse translate reads into three .msa files. This function takes as arguments
# a path to a directory containing the input files and a path to a directory in
# which the output files will be written.

def collapse_translate_reads (in_path, out_path):
    
    # Identify the paths to the three output files.
    
    merged_out = join (out_path, "merged.msa")
    translated_prot = join (out_path, "traslated.msa")
    merged_out_prot = join (out_path, "merged_prot.msa")
    
    # Check that these paths exist.
    
    if check_file_paths_diversity and os.path.exists (merged_out) and  os.path.exists (merged_out_prot):
        return (merged_out, merged_out_prot)  
    
    print ("Collapsing and translating reads %s " % in_path, file = sys.stderr)
    
    # Call the seqmerge and translate utilities on the input files, then call
    # seqmerge on the translate utility's output.
    
    try:
        subprocess.check_call (['/opt/share/python3.3/seqmerge', in_path, merged_out], stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL) 
        subprocess.check_call (['/opt/share/python3.3/translate', in_path, translated_prot], stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL) 
        subprocess.check_call (['/opt/share/python3.3/seqmerge', translated_prot, merged_out_prot], stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
        os.remove (translated_prot)
    
    # If an error occurs, say so.
    
    except subprocess.CalledProcessError as e:
        print ('ERROR: Collapse/translate call failed failed',e,file = sys.stderr)
        return None
        
    # Return the paths to the merged output files.
    
    return (merged_out, merged_out_prot)   




### count_collapsed_reads ###

# Count the collapsed read files. This function takes as agruments a path to
# a directory containing the collapsed files, a path to a directory in which
# a prot coverage .json file will be recorded, and a node to run on.

def count_collapsed_reads (in_path, out_path, node):
    print ("Getting protein coverage info for %s (node %d) " % (in_path, node), file = sys.stderr)
    
    # Identify the output file.
    
    merged_json= join (out_path, "prot_coverage.json")
    
    # Check that this file path exists, and load any .json cache that has
    # already been written there.

    if check_file_paths_diversity and os.path.exists (merged_json):
        with open (merged_json) as fh:
            if len(json.load(fh)) > 0: 
                print ("[CACHED]")
                return merged_json 
                
    # Use bpsh to call the seqcoverage utility on the coverage .json file,
    # using the selected node. 
    
    try:
        #print (' '.join(['/usr/bin/bpsh', str(node), '/usr/local/bin/seqcoverage', '-o', merged_json, '-t', 'protein', in_path]))
        subprocess.check_call (['/usr/bin/bpsh', str(node), '/usr/local/bin/seqcoverage', '-o', merged_json, '-t', 'protein', in_path], stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL) 
        
    # If an error occurs, say so.
    
    except subprocess.CalledProcessError as e:
        print ('ERROR: Protein coverage call failed',e,file = sys.stderr)
        return None
        
    # Return the path to the .json cache file.
    
    return merged_json   




### multinomial_filter ###

# apply multinomial filter

def multinomial_filter (in_path, out_path, node):
    filtered_out = join (out_path, "filtered.msa")
    json_out = join (out_path, "rates.json")

    if check_file_paths and os.path.exists (filtered_out) and os.path.exists (json_out):
        return filtered_out, json_out
    try:
        print ("Running multinomial filter on %s (node = %d) " % (in_path, node), file = sys.stderr)
        #-t 0.005 -p 0.999999 -f results/Ionxpress020/filtered.msa -j results/Ionxpress020/rates.json results/Ionxpress020/aligned.msa
        subprocess.check_call (['/usr/bin/bpsh', str(node), os.path.join (path_to_this_file, "../julia/mcmc.jl"), '-t', '0.005', '-p', '0.999', '-f', filtered_out, '-j', json_out, in_path], stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL) 
    except subprocess.CalledProcessError as e:
        print ('ERROR: multinomial filter call failed',e,file = sys.stderr)
        return None, None
    return filtered_out, json_out    

# check compartmentalization using tn93
    
def check_compartmenalization (in_paths, node, replicates = 100, subset = 0.2, min_overlap = 150):
    print ("Running compartmenalization tests on %s (node %d) " % (in_paths, node), file = sys.stderr)
    baseline_json = None
    out = ''
    json_out = ''
    try:
        status = 'Running initial F_ST'
        process = subprocess.Popen (['/usr/bin/bpsh', str(node), '/usr/local/bin/tn93', '-t', str(0.01), '-l', str (min_overlap), '-c', '-q', '-m', '-u', str(subset), '-s', in_paths[0][0], in_paths[1][0]], stdin = subprocess.DEVNULL, stderr = subprocess.PIPE, stdout = subprocess.PIPE, universal_newlines = True) 
        out, json_out = process.communicate ()
        baseline_json = json.loads(json_out)
        baseline = baseline_json ['F_ST']
        #print ("F_ST baseline = %g" % baseline)
        p_v = 0
        for k in range (replicates):
            status = 'Running replicate %d' % k
            process = subprocess.Popen (['/usr/bin/bpsh', str(node), '/usr/local/bin/tn93', '-t', str(0.01), '-l', str (min_overlap), '-c', '-b', '-q', '-m', '-u', str(subset), '-s', in_paths[0][0], in_paths[1][0]], stdin = subprocess.DEVNULL, stderr = subprocess.PIPE, stdout = subprocess.PIPE, universal_newlines = True) 
            out, json_out = process.communicate ()
            sim_fst = json.loads(json_out)['F_ST']
            p_v += 1 if sim_fst >= baseline else 0
            #print ("%d %g %g" % (k, p_v/replicates, sim_fst))
        
        p_v = (p_v+1.)/(replicates+1.)
            
        baseline_json = {in_paths[0][1] : baseline_json['Histogram File 1'], 
                         in_paths[1][1] : baseline_json['Histogram File 2'], 
                         'Between': baseline_json['Histogram Between'], 
                         'p' : p_v,
                         'f_st' : baseline}
                    
        return baseline_json
                
    except subprocess.CalledProcessError as e:
        print ('ERROR: tn93 call failed in check_compartmenalization',e,file = sys.stderr)
    except Exception as e:
        print ('ERROR: check_compartmenalization error on %s and %s\n\t\n%s\n%s\n%s' % (str(in_paths[0]), str(in_paths[1]), status, json_out, out),e,file = sys.stderr)
        pass
        #print (json_out)

    return baseline_json

# extract diagnostic region 

def extract_diagnostic_region (in_path, out_path, node, start = 0, end = 1000000, cov = 0.95):
    #print ("Extracting diagnostic region [%d-%d, coverage = %g] for %s (node %d) " % (start, end, cov, in_path, node), file = sys.stderr)
    merged_out = join (out_path, "region_%d-%d.msa" % (start, end))
    
    if check_file_paths_diversity and os.path.exists (merged_out) and not force_diversity_estimation:
        return merged_out
    try:
        call_array = ['/usr/bin/bpsh', str(node), '/usr/local/bin/selectreads', '-o', merged_out, '-a', 'gaponly', '-s', str (start), '-e', str (end), '-c', str (cov), in_path]
        #print (call_array)
        subprocess.check_call (call_array,stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL) 
    except subprocess.CalledProcessError as e:
        print ('ERROR: Diagnostic region extract call failed',e,file = sys.stderr)
        return None
    return merged_out      
    
# collapse_diagnostic_region

def collapse_diagnostic_region (in_path_list, out_path, node, overlap = 100, count = 32):
    result = {}
    for region, in_path in in_path_list.items():    
        if in_path is not None and os.path.exists (in_path) or force_diversity_estimation:
            merged_out = join (out_path, "region_reduced_%s.msa" % region)
            try:
                if check_file_paths_diversity and os.path.exists (merged_out) and not force_diversity_estimation:
                    result [region] = merged_out
                    continue
                print ("Collapsing diagnostic region (overlap %d, count = %d) for %s (node %d) " % (overlap, count, in_path, node), file = sys.stderr)
                call_array = ['/usr/bin/bpsh', str(node), '/usr/local/bin/readreduce', '-o', merged_out, '-l', str (overlap), '-s', str (count), in_path]
                #print (call_array)
                subprocess.check_call (call_array, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL) 
                #subprocess.check_call (call_array, stdout = subprocess.DEVNULL) 
                if os.stat (merged_out).st_size > 0:
                    result [region] = merged_out
                else:
                    os.remove (merged_out)
                
                os.remove (in_path)
            
            except subprocess.CalledProcessError as e:
                print ('ERROR: Diagnostic region collapse call failed',e,file = sys.stderr)
                return None
            except:
                raise
                
    return result         

# load a merged file and find, extract, and collapse well covered regions

def extract_and_collapse_well_covered_region (in_path, out_path, node, read_length = 200, min_coverage = 100, min_read_count = 10):
    merged_out = join (out_path, "prot_coverage.json")
    with open (merged_out) as fh:
        coverage = json.load (fh)
        positional_coverage = [[int (k), sum (v.values())] for k, v in coverage.items()]
        positional_coverage = [k for k in positional_coverage if k[1] >= min_coverage]
        if len (positional_coverage) == 0:
            return None
        median = max(min_coverage, describe_vector ([ k [1] for k in positional_coverage])['median'] * 0.5)
        sorted_positions = sorted([ k for k in positional_coverage if k[1] >= median])
        print ("Extracting well-covered region for (read-length = %g, min_coverage = %d, selected_coverage = %g) %s (node %d) " % (float(read_length), min_coverage, median, in_path, node), file = sys.stderr)
        if len (sorted_positions) > 50:
            '''
            end   = len(sorted_positions) // 2
            start = end
            while start > 0  and sorted_positions[start-1][0] == sorted_positions[start][0] - 1:
                start -= 1
            while end < len (sorted_positions)-1 and sorted_positions[end+1][0] == sorted_positions[end][0] + 1:
                end += 1
            '''
                          
            start = sorted_positions[0][0] - 1
            end = sorted_positions[-1][0] - 1
            print ("Selected [%d - %d]" % (start, end), file = sys.stderr)
            #print (sorted_positions)
            
            if end - start >= 50:
                start = start*3
                end   = end*3
                region_key = str(start) + "-" + str(end)
                extracted = extract_diagnostic_region (in_path, out_path, node, start, end, cov = min (0.95, float(max(window,read_length))/(end-start)))
                collapsed = collapse_diagnostic_region ({region_key:extracted}, out_path, node, overlap = max (50, float (read_length * 0.50)), count = max(10, int (median * 0.005)))   
                if collapsed is not None:
                    if region_key in collapsed:       
                        return list(collapsed.items())[0]
                
    return None

# run tropism prediction

def run_tropism_prediction (env_gene, out_path, node):
    result = {}
    v3_model = os.path.join (path_to_this_file, "../data/V3.model")
    
    idepi_out = join (out_path, "idepi.json")
    updated = False
    
    if os.path.exists (idepi_out):
        with open (idepi_out, "r") as fh:
            try:
                result = json.load (fh)
            except ValueError as e:
                pass
                #print ("Error reloading JSON info from %s" % diversity_out, file = sys.stderr)
                #raise e
                
    if len (result) == 0:
        print ("Running tropism predictions on %s (node %d)" % (env_gene, node), file = sys.stderr)
        with open (idepi_out, "w") as json_out:
            try:
                call_array = ['/usr/bin/bpsh', str(node), '/opt/share/python3.3/idepi', 'predict', v3_model, env_gene]
            #print (call_array)
                subprocess.check_call (call_array,stdout = json_out, stderr = subprocess.DEVNULL) 
            except subprocess.CalledProcessError as e:
                print ('ERROR: Tropism predictions call failed',e,file = sys.stderr)
                return None

        with open (idepi_out, "r") as fh:
            try:
                result = json.load (fh)
            except ValueError as e:
                pass
        
    if len (result):
        counts = {-1 : 0, 1 : 0}
        
        for seq in result['predictions']:
            counts [seq['value']] += float(seq['id'].split (':')[1])
            
        total = sum (counts.values())
        return {'R5' :  counts[-1] / total, 'X4' : counts[1] / total}
    
    return None

# process an extracted diagnostic region

def process_diagnostic_region (in_path_list, out_path, node):
    result = {}
    file = os.path.join (path_to_this_file, "../hyphy", "maketree.bf")
    
    diversity_out = join (out_path, "diversity.json")
    updated = False
    
    if os.path.exists (diversity_out) and not force_diversity_estimation:
        with open (diversity_out, "r") as fh:
            try:
                result = json.load (fh)
            except ValueError as e:
                #print ("Error reloading JSON info from %s" % diversity_out, file = sys.stderr)
                #raise e
                result = {}

    for region, in_path in in_path_list.items():    
        if region in result:
            continue
            
        print ("Running diversity estimates on %s (node %d)" % (in_path, node), file = sys.stderr)
        try:
            in_file = os.path.join (path_to_this_file, in_path)
            process = subprocess.Popen (['/usr/bin/bpsh', str(node), '/usr/local/bin/HYPHYMP', file], stdin = subprocess.PIPE, stderr = subprocess.PIPE, stdout = subprocess.PIPE) 
            out, err = process.communicate (bytes(in_file,'UTF-8'))
            try:
                out = json.loads(out.decode ('UTF-8'))
                result [region] = out
                updated = True
            except ValueError:
                print ('ERROR HyPhy call (file %s) %s/%s' % (in_file, out.decode ('UTF-8'), err.decode ('UTF-8')), file = sys.stderr)
                pass
            
        except subprocess.CalledProcessError as e:
            print ('ERROR: Diagnostic region analysis call failed',e,file = sys.stderr)
            return None, False
            
    with open (diversity_out, "w") as fh:
        json.dump (result, fh, sort_keys=True, indent=1)
               
    return diversity_out, updated    
    
## END TASK RUNNERS ##      

def update_global_record (base_path, gene, analysis_record):
    global threading_lock
    global NGS_run_cache

    threading_lock.acquire()
    NGS_run_cache[base_path][gene] = analysis_record
    print ("DUMPING JSON WITH %d RECORDS" % len (NGS_run_cache))
    with open (cache_file, "w") as fh:
        json.dump (NGS_run_cache, fh,default=dthandler, sort_keys=True, indent=1)
    update_json = False
    
    threading_lock.release()
        
def analysis_handler (node_to_run_on):
    while True:
        base_path, file_results_dir_overall, i, gene, analysis_cache, median_read_length = task_queue.get()
        handle_a_gene (base_path, file_results_dir_overall, i, gene, analysis_cache, node_to_run_on, median_read_length)
        task_queue.task_done()

def compartmentalization_handler (node_to_run_on):
    while True:
        in_paths, tag, analysis_record = task_queue.get()
        cr = check_compartmenalization (in_paths, node_to_run_on, subset = 0.33)
        threading_lock.acquire()
        analysis_record[tag] = cr
        
        with open (cache_file, "w") as fh:
            json.dump (NGS_run_cache, fh,default=dthandler, sort_keys=True, indent=1)
            
        threading_lock.release()
        task_queue.task_done()
 
 
def set_update_json (path, text):
    print ("Updated JSON for %s because of %s" % (path, text), file = sys.stderr)
    return True

def handle_a_gene (base_path, file_results_dir_overall, i, gene, analysis_cache, node, median_read_length):
     
    global threading_lock
    global NGS_run_cache
    
    update_json = False
    file_results_dir = os.path.join (file_results_dir_overall, gene)

    if not os.path.exists (file_results_dir):
        os.makedirs (file_results_dir) 
        
    if i == 0:
        if 'filtered_fastq' in NGS_run_cache [base_path] and NGS_run_cache [base_path] ['filtered_fastq'] is not None:
            if 'aligned_bam' not in analysis_cache or analysis_cache['aligned_bam'] is None:
                analysis_cache ['aligned_bam'], analysis_cache ['discards'] = ntr.codon_aligner (NGS_run_cache [base_path] ['filtered_fastq'], file_results_dir, node, refs[i])
                update_json = True
            if 'reference_sequence' not in analysis_cache:
                analysis_cache['reference_sequence'] = refs[i]
                update_json = set_update_json (file_results_dir, "if 'reference_sequence' not in analysis_cache")
    else:
        prev_gene = genes[i-1]
        if 'discards' in NGS_run_cache [base_path][prev_gene] and NGS_run_cache [base_path][prev_gene]['discards'] is not None:
            if 'aligned_bam' not in analysis_cache or analysis_cache['aligned_bam'] is None:
                analysis_cache ['aligned_bam'], analysis_cache ['discards'] = ntr.codon_aligner (NGS_run_cache [base_path][prev_gene]['discards'], file_results_dir, node, refs[i])
                update_json = True
            if 'reference_sequence' not in analysis_cache:
                analysis_cache['reference_sequence'] = refs[i]
                update_json = set_update_json (file_results_dir, "if 'reference_sequence' not in analysis_cache")



    if update_json:
        update_global_record (base_path, gene, analysis_cache)
        update_json = False

    if 'aligned_bam' in analysis_cache and analysis_cache ['aligned_bam'] is not None:
        if 'aligned_msa' not in analysis_cache or analysis_cache['aligned_msa'] is None:
            analysis_cache ['aligned_msa'] = ntr.bam_to_fasta (analysis_cache ['aligned_bam'], file_results_dir)
            if analysis_cache ['aligned_msa'] is not None and os.path.getsize ( analysis_cache ['aligned_msa']) > 0:
                update_json = set_update_json (file_results_dir, "os.path.getsize ( analysis_cache ['aligned_msa']) > 0")
            else:
                analysis_cache['aligned_msa'] = None

    if 'aligned_msa' in analysis_cache and analysis_cache ['aligned_msa'] is not None:
        if 'filtered_msa' not in analysis_cache or 'json_rates' not in analysis_cache:
            analysis_cache ['filtered_msa'], analysis_cache ['json_rates']  = multinomial_filter (analysis_cache ['aligned_msa'], file_results_dir, node)
            update_json = set_update_json (file_results_dir, "'filtered_msa' not in analysis_cache or 'json_rates' not in analysis_cache")
            
 
    if update_json:
        update_global_record (base_path, gene, analysis_cache)
        update_json = False

    if 'filtered_msa' in analysis_cache and analysis_cache ['filtered_msa'] is not None:
        if 'merged_msa' not in analysis_cache or 'merged_msa_prot' not in analysis_cache:
            analysis_cache ['merged_msa'], analysis_cache ['merged_msa_prot']  = collapse_translate_reads (analysis_cache ['filtered_msa'], file_results_dir)
        if force_diversity_estimation or 'region_msa' not in analysis_cache:
            result = {}
            for span in spans[i]:
                result["-".join([str (k) for k in span])] = extract_diagnostic_region (analysis_cache ['merged_msa'], file_results_dir, node, span[0], span[1])
            
            
            analysis_cache['region_msa'] = result
            update_json = set_update_json (file_results_dir, "'region_msa' not in analysis_cache")

    if 'merged_msa_prot' in analysis_cache and analysis_cache ['merged_msa_prot'] is not None:
        if 'merged_json' not in analysis_cache:
            analysis_cache ['merged_json'] = count_collapsed_reads (analysis_cache ['merged_msa_prot'], file_results_dir, node)
            update_json = set_update_json (base_path, "'merged_json' not in analysis_cache")
        else:
            with open (analysis_cache ['merged_json'], 'r') as fh:
                if len(json.load (fh)) == 0:
                    analysis_cache ['merged_json'] = count_collapsed_reads (analysis_cache ['merged_msa_prot'], file_results_dir, node)
                    update_json = set_update_json (file_results_dir, "len(json.load (fh)) == 0")
                    
            
        if force_diversity_estimation or 'overall_region' not in analysis_cache or analysis_cache['overall_region'] is None:
            analysis_cache['overall_region'] = extract_and_collapse_well_covered_region (analysis_cache ['merged_msa'], file_results_dir, node, median_read_length) 
            update_json = set_update_json (file_results_dir, "'overall_region' not in analysis_cache or analysis_cache['overall_region'] is None")
                                                      

    if update_json:
        update_global_record (base_path, gene, analysis_cache)
        update_json = False
        
    if 'region_msa' in analysis_cache and analysis_cache ['region_msa']  is not None:
        if force_diversity_estimation or 'region_merged_msa' not in NGS_run_cache [base_path] [gene] or NGS_run_cache [base_path] [gene]['region_merged_msa'] is None:
            analysis_cache ['region_merged_msa'] = collapse_diagnostic_region (analysis_cache ['region_msa'], file_results_dir, node, window - 10)
            update_json = set_update_json (file_results_dir, "'region_merged_msa' not in NGS_run_cache [base_path] [gene] or NGS_run_cache [base_path] [gene]['region_merged_msa'] is None")

    if update_json:
        update_global_record (base_path, gene, analysis_cache)
        update_json = False

    if 'region_merged_msa' in analysis_cache and analysis_cache ['region_merged_msa']   is not None:
        if 'overall_region' in analysis_cache and analysis_cache ['overall_region'] is not None:
            if force_diversity_estimation or analysis_cache ['overall_region'][0] not in analysis_cache['region_merged_msa']:
                analysis_cache['region_merged_msa'][analysis_cache ['overall_region'][0]] = analysis_cache ['overall_region'][1]
                update_json = set_update_json (file_results_dir, "analysis_cache ['overall_region'][0] not in analysis_cache['region_msa']")
            
            if gene == 'env' and ('tropism' not in analysis_cache or analysis_cache['tropism'] is None):
                analysis_cache['tropism']  = run_tropism_prediction (analysis_cache ['overall_region'][1], file_results_dir, node)
                
        analysis_cache ['region_merged_processed'], update_json_local = process_diagnostic_region (analysis_cache ['region_merged_msa'], file_results_dir, node)
        if update_json_local:
            update_json = set_update_json (file_results_dir, "process_diagnostic_region (analysis_cache ['region_merged_msa'], file_results_dir, node)")

    if 'merged_msa' in analysis_cache and analysis_cache ['merged_msa'] is not None:
        if 'tn93_json' not in analysis_cache:
            analysis_cache ['tn93_json'] = ntr.get_tn93 (analysis_cache ['merged_msa'], file_results_dir, node)
            update_json = set_update_json (file_results_dir, "'merged_msa' in analysis_cache and analysis_cache ['merged_msa'] is not None")
        
        '''
        if has_compartment_data:
            tag = (patient_id, sample_date)
            if tag not in check_for_compartmentalization:
                check_for_compartmentalization[tag] = {'genes' : [], 'paths' : {}}
        
            if gene not in check_for_compartmentalization[tag]['genes']:
                check_for_compartmentalization[tag]['genes'].append (gene)
                check_for_compartmentalization[tag]['paths'][gene] = []
            
            check_for_compartmentalization[tag]['paths'][gene].append ([file_results_dir,analysis_cache['merged_msa'],compartment,base_path])    
        '''
    if update_json:
        update_global_record (base_path, gene, analysis_cache)
        
def describe_vector (vector):
    vector.sort()
    l = len (vector)
    if l:
        return {'count': l, 'min': vector[0], 'max': vector[-1], 'mean': sum(vector)/l, 'median':  vector [l//2] if l % 2 == 1 else 0.5*(vector[l//2-1]+vector[l//2]), "IQR": [vector [l//4], vector [(3*l)//4]] }
    else:
        return {'count': l, 'min': None, 'max': None, 'mean': None, 'median':  None, "IQR": [None,None] }


def check_keys_in_dict (d, k):
    if k[0] not in d:
        d [k[0]] = {}
    if len (k) > 1:
        return check_keys_in_dict (d [k[0]], k[1:])   
    return d[k[0]]
    
def hash_file (filepath):
    m = hashlib.md5()
    with open (filepath, "rb") as fh:
        m.update (fh.read())
    return m.hexdigest() 
    
def main (dir, results_dir, has_compartment_data, has_replicate_counts, scan_q_filt):
                
    global NGS_run_cache
    global task_queue
    
    non_gene_keys = set (('id','in_fasta','in_qual','patient_id','compartment','filtered_fastq','replicate','sample_date'))                 

    if os.path.exists(cache_file):
        with open (cache_file, "r") as fh:
            NGS_run_cache = json.load (fh)
    else:
        NGS_run_cache = {}
        

            
    #refs  = ['data/pr.fas', 'data/rt.fas','data/gag_p24.fas','data/env_C2V5.fas']
    #genes = ['pr','rt','gag', 'env']
  
    #for k in spans:            
    #    print (k)

    
    tasks_by_gene = {}

    # check for zip files first
    for root, dirs, files in os.walk(dir):
        for file in files:
            name, ext = splitext (file)
            if len (ext) > 0 and ext in ('.zip'):
                base_path = join(root, name) 
                full_file = base_path + ".zip" 
                print ("Unzipping %s " % full_file)
                try:
                    subprocess.check_call (['/usr/bin/unzip', '-n', '-d', root, '-j', full_file])
                    os.remove (full_file)            
                except subprocess.CalledProcessError as e:
                    print ('ERROR: UNZIP call failed failed',e,file = sys.stderr)
                
    #sys.exit (0)
        
    
    for root, dirs, files in os.walk(dir):
        for file in files:
            name, ext = splitext (file)
            if len (ext) > 0 and ext in ('.fna','.fastq'):
                base_path = root  
                if name == 'qfilt' and scan_q_filt != True or name == 'discards': 
                    continue
                base_file = join(root, name + ext)

                #threading_lock.acquire()
                if base_path not in NGS_run_cache:
                    NGS_run_cache [base_path] = {'id' : len (NGS_run_cache) + 1}
                    
                print ('Working on %s...' % base_file, file = sys.stderr)             
                
                
                if has_replicate_counts:
                    if has_compartment_data:
                        dirn, replicate = os.path.split (root)    
                        dirn, compartment = os.path.split (dirn)    
                        dirn, sample_date = os.path.split (dirn)
                        dirn, patient_id  = os.path.split (dirn)
                        file_results_dir_overall = os.path.join (results_dir, patient_id, sample_date, compartment, replicate)
                        NGS_run_cache[base_path]['patient_id'] = patient_id
                        NGS_run_cache[base_path]['sample_date'] = sample_date
                        NGS_run_cache[base_path]['compartment'] = compartment
                        NGS_run_cache[base_path]['replicate'] = replicate
                    else:            
                        dirn, replicate = os.path.split (root)    
                        dirn, sample_date = os.path.split (dirn)
                        dirn, patient_id = os.path.split (dirn)
                        file_results_dir_overall = os.path.join (results_dir, patient_id, sample_date, replicate)                
                        NGS_run_cache[base_path]['patient_id'] = patient_id
                        NGS_run_cache[base_path]['sample_date'] = sample_date
                        NGS_run_cache[base_path]['replicate'] = replicate
                else:
                    if has_compartment_data:
                        dirn, compartment = os.path.split (root)    
                        dirn, sample_date = os.path.split (dirn)
                        dirn, patient_id  = os.path.split (dirn)
                        file_results_dir_overall = os.path.join (results_dir, patient_id, sample_date, compartment)
                        NGS_run_cache[base_path]['patient_id'] = patient_id
                        NGS_run_cache[base_path]['sample_date'] = sample_date
                        NGS_run_cache[base_path]['compartment'] = compartment
                    else:            
                        dirn, sample_date = os.path.split (root)
                        dirn, patient_id = os.path.split (dirn)
                        file_results_dir_overall = os.path.join (results_dir, patient_id, sample_date)
                        NGS_run_cache[base_path]['patient_id'] = patient_id
                        NGS_run_cache[base_path]['sample_date'] = sample_date
                                
                                
                if not os.path.exists (file_results_dir_overall):
                    os.makedirs (file_results_dir_overall) 

                if 'in_fasta' not in NGS_run_cache [base_path] and 'in_fastq' not in NGS_run_cache [base_path]:
                    NGS_run_cache [base_path]['md5'] = hash_file (base_file)                                         
                    if ext == '.fastq':
                        NGS_run_cache [base_path]['in_fastq'] = base_file  
                    else:
                        NGS_run_cache [base_path]['in_fasta'] = base_file
                        NGS_run_cache [base_path]['in_qual'] = join(root, name + ".qual")
                else:
                    do_skip  = False
                    same_hash = False
                    
                    key_pair = ('in_fastq', '.fna') if 'in_fastq' in NGS_run_cache [base_path] else ('in_fasta', '.fastq')
                    
                    if ext == key_pair[1]:
                        do_skip = True
                    else:
                        if NGS_run_cache [base_path][key_pair[0]] != base_file:
                            do_skip = True
                            same_hash =  NGS_run_cache [base_path]['md5'] == hash_file (base_file)
                    
                    if do_skip:
                        print ("Skipping file %s because the path has already been processed, i.e. there are multiple NGS files in %s'" % (base_file, root) , file = sys.stderr)
                        if not same_hash:
                            print ("\tWARNING!!!! Different hashes for different NGS files")
                        continue
    

                median_read_length = 200
                
                if 'in_fasta' in NGS_run_cache [base_path] and  NGS_run_cache [base_path] ['in_fasta'] is not None:
                    if 'filtered_fastq' not in NGS_run_cache [base_path] or force_qfilt_rerun:
                        NGS_run_cache [base_path] ['filtered_fastq'] = run_qfilt (NGS_run_cache [base_path] ['in_fasta'], NGS_run_cache [base_path] ['in_qual'], join (file_results_dir_overall, 'qfilt.fna'), join (file_results_dir_overall, 'qfilt.json'))

                if 'in_fastq' in NGS_run_cache [base_path] and  NGS_run_cache [base_path] ['in_fastq'] is not None:
                    if 'filtered_fastq' not in NGS_run_cache [base_path] or force_qfilt_rerun:
                        NGS_run_cache [base_path] ['filtered_fastq'] = run_qfilt (NGS_run_cache [base_path] ['in_fastq'], None, join (file_results_dir_overall, 'qfilt.fna'), join (file_results_dir_overall, 'qfilt.json'))
                    
                if 'filtered_fastq' in NGS_run_cache [base_path] and NGS_run_cache [base_path]['filtered_fastq'] is not None:
                    try:
                        with open (join (file_results_dir_overall, 'qfilt.json')) as fh: 
                            js = json.load (fh)['run summary']['original read length distribution:']
                            median_read_length = max (100, js['mean'] - js['standard deviation'])
                    except:
                        pass
                        
                #threading_lock.release()

                for i, gene in enumerate(genes): 
                    if gene not in NGS_run_cache [base_path]:
                        NGS_run_cache [base_path][gene] = {} 
                    if gene not in tasks_by_gene:
                        tasks_by_gene[gene] = []
                    tasks_by_gene[gene].append ([base_path, file_results_dir_overall, i, gene, NGS_run_cache [base_path][gene], median_read_length])
                    

    
    with open (cache_file, "w") as fh:
        json.dump (NGS_run_cache, fh,default=dthandler, sort_keys=True, indent=1)
    
    for gene,task_list in tasks_by_gene.items():
        task_queue = queue.Queue ()
    
        for node in nodes_to_run_on:
            t = threading.Thread (target = analysis_handler, args = (node,))
            t.daemon = True
            t.start ()

        for task in task_list:
            task_queue.put (task, block = False)
            
        task_queue.join()    
    #task_queue.join()    

    task_queue = queue.Queue ()
    if 'F_ST' not in NGS_run_cache:
        NGS_run_cache['F_ST'] = {}
        
    compartmentalization_sets = {}
    
    for key,value in NGS_run_cache.items():
        if type (value) == dict:
            if 'patient_id' in value and 'sample_date' in value and 'compartment' in value:
                tag = (value['patient_id'], value['sample_date'])
                if tag not in compartmentalization_sets:
                    compartmentalization_sets[tag] = {}
                for key2, value2 in value.items():
                    if key2 not in non_gene_keys:
                        if 'merged_msa' in value2 and 'overall_region' in value2 and value2['overall_region']:
                            #print (tag, key2, value['compartment'] )
                            if key2 not in compartmentalization_sets[tag]:
                                compartmentalization_sets[tag][key2] = {}
                            if value['compartment'] not in compartmentalization_sets[tag][key2]:
                                compartmentalization_sets[tag][key2][value['compartment']] = []
                            compartmentalization_sets[tag][key2][value['compartment']].append ( value2['merged_msa'])
    
    for node in nodes_to_run_on:
        t = threading.Thread (target = compartmentalization_handler, args = (node,))
        t.daemon = True
        t.start ()

    for sample, info in compartmentalization_sets.items():
        for gene, data in info.items():
            if len (data) >= 2:
                compartments = list (data.keys())
                
                subject_cache = check_keys_in_dict (NGS_run_cache['F_ST'], [sample[0], sample[1], gene])
                                
                                
                for comp1, comp2 in combinations (compartments, 2):
                    for f1, f2 in product (data[comp1], data[comp2]):
                        pair_tag = "%s|%s" % ((f1, f2) if f1 < f2  else (f2, f1))
                        if pair_tag not in subject_cache or subject_cache[pair_tag] is None or 'p' not in subject_cache[pair_tag]:
                            task_queue.put ([[[f1, comp1], [f2, comp2]],pair_tag,subject_cache])
                            #subject_cache [pair_tag] = check_compartmenalization ([[f1, compartments[0]], [f2, compartments[1]]], 2, subset = 0.2)
                            #print (subject_cache[pair_tag])
                            #task_queue.join()    
                            #sys.exit (1)
                        #else:
                        #    print (subject_cache[pair_tag])
                        
                
                    
    task_queue.join()    

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
    
    parser.add_argument(
        '-p', '--compartment',
        action='store_true',
        help='does the directory structure include compartment information',
    )    

    parser.add_argument(
        '-d', '--replicate',
        action='store_true',
        help='does the directory structure include replicate information',
    )    

    parser.add_argument(
        '-q', '--qfilt',
        action='store_true',
        help='scan qfilt.fna files',
    )    

    parser.add_argument(
        '-n', '--node',
        type=str,
        required=True,
        help='run MP scripts on this node',
        default = "1"
    )    

    parser.add_argument(
            '-g', '--genes',
            metavar='genes',
            type=str,
            help='the comma separated list of genes to include in the comparison',
            default='env,gag,rt'
    )

    
    threading_lock = threading.Lock()
    
    args = None
    retcode = -1
    args = parser.parse_args()
    
    if not os.path.exists (args.results):
        os.mkdir (args.results) 
    
    nodes_to_run_on = [int (k) for k in args.node.split (",")]
    
    genes = args.genes.split (',')
    
    try:
        refs = [known_refs [known_genes.index (k)] for k in genes] 
    except:
        print ('Please check that all the genes come from the following list %s' % str (known_genes))
        sys.exit (1)
    
    
    for i,file in enumerate(refs):
        with open (file) as fh:
            for record in SeqIO.parse(fh, "fasta") :
                reflen = len(record.seq)
                spans.append([[k,k+window] for k in range (0,reflen-window,stride)])
       
    #force_diversity_estimation = True
    #extract_and_collapse_well_covered_region ("/data/collaborators/gert_van_zyl/A144/20080513/miSeqR2/rt/merged.msa", "/data/collaborators/gert_van_zyl/A144/20080513/miSeqR2/rt", 2)
    #sys.exit (retcode)
    
    cache_file = args.cache
    retcode = main(args.input,  args.results, args.compartment, args.replicate, args.qfilt)
    
    sys.exit(retcode)
