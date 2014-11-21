import os, argparse, re, sys, datetime
import subprocess, time, csv
import datetime, shutil, json, operator
from os.path import join, splitext, isfile, split
from Bio import SeqIO, AlignIO
from copy import copy
from scipy import stats
import numpy as np
import math

DRAMs = {}

def parse_a_date (sample_date):
    try:
        return datetime.datetime.strptime(sample_date, "%Y%m%d").strftime ("%Y/%m/%d")
    except ValueError as e:
        return sample_date


def add_key_if_missing (d, k, v):
    if k not in d:
        d[k] = v
        
    return d [k]

def describe_vector (vector):
    vector.sort()
    l = len (vector)
    return {'count': l, 'min': vector[0], 'max': vector[-1], 'mean': sum(vector)/l, 'median':  vector [l//2] if l % 2 == 1 else 0.5*(vector[l//2-1]+vector[l//2]), "IQR": [vector [l//4], vector [(3*l)//4]] }


def coverage_info (coverage, min_cov = 500):
    ranges = sorted([int (k) for k in coverage])

    if (len (coverage) == 0):
        return (None, None)

    region_coverage = []
    unfiltered    = []
    
    
    for s in range (ranges[0], ranges[-1]+1):
        i = str (s)
        if i in coverage:
            val = sum (coverage[i].values())
            unfiltered.append (val)
            if val >= min_cov:
                region_coverage.append ([s,val])
            
        
    coordinates = [k[0] for k in region_coverage]
    coverages   = [k[1] for k in region_coverage]
    
    return  ([min(coordinates), max(coordinates)], describe_vector(coverages)) if len(coordinates) > 0 else (None, describe_vector(unfiltered))

def diversity_info (div, tag):
    try:
        values = [div[k]["Diversity"][tag] for k in div]
        return    describe_vector(values)
    except:
        return None


def main (cache_file, out, store_dir, has_compartment, has_replicate, intrahost_info, pairwise_distance_info, short, ignore_replicate):
    
    has_tropism = False
    
    processed    = {'columns': ['PID','Date','Gene','Region','Median Coverage','Nucleotide diversity, %','Syn.diversity, %', 'Non-syn. diversity, %' ,'X4 tropism, %',''], 
                      'data' : [] }
                      
    if intrahost_info is not None:
        processed ['intrahost'] = {}
        
    has_fst = "F_ST" in cache_file

    processed['settings'] = {'compartment' : has_compartment, 'replicate': has_replicate}
    
    if has_replicate:
        processed['columns'].insert (3, 'Replicate')
    if has_compartment:
        processed['columns'].insert (3, 'Compartment')
        
    if has_fst:
        processed['F_ST'] = {}
        
                                    
    '''
    for root, dirs, files in os.walk(dir):
        for d in dirs:
            available_information [os.path.join (root, d)] = {'prot_coverage' : None, 'diversity' : None, 'tn93' : None, 'rates': None}
            optional_information [os.path.join (root, d)]  = {'fst' : None}
            
        for file in files:
            name, ext = splitext (file)
            for which_dict in (available_information, optional_information):
                if root in which_dict and name in which_dict[root]:
                    try:
                        with open (join (root, file)) as fh:
                            which_dict[root][name] = json.load (fh)
                    except:
                        pass
    '''                
 
    required_keys = ['merged_json', 'region_merged_processed', 'tn93_json', 'json_rates']
    optional_keys = ['fst']
    
    id_to_pid_date       = {}
    earilest_date_by_pid = {}
    merged_msa_to_info   = {}
    
    print ("Loaded the pipeline cache file with %d records" % len (cache_file), file = sys.stderr)
    
    tried = 0
    
    for record, item in cache_file.items():
        try:
            if 'filtered_fastq' not in item or item['filtered_fastq'] is None:
                continue

            row = []
            pid  = item ["patient_id"]
            sample_date = item["sample_date"]
            date = parse_a_date (sample_date)
            
                
            if pid not in earilest_date_by_pid:
                earilest_date_by_pid[pid] = date
            else:
                earilest_date_by_pid[pid] = min (date, earilest_date_by_pid[pid])
            
            id_expansion = [pid, date]
            
            if has_compartment:
                compartment = item["compartment"]
                id_expansion.append (compartment)
            else:
                compartment = None
                
            if has_replicate:
                replicate = item["replicate"]
                id_expansion.append (replicate)
            else:
                replicate = None
            id_to_pid_date [item["id"]] = id_expansion
            
            
                
            for gene, data in item.items():                     
                try:
                    current_step = None
                    dir_name = [pid, sample_date, gene]
                    
                    tried += 1
                
                    row = [pid,date,gene]
                    if has_compartment:
                        row.insert (3,compartment )
                        dir_name.append (compartment)
                        
                    if has_replicate:
                        row.insert (4, replicate)
                        dir_name.append (replicate)
                      
                    dir_name = '_'.join (dir_name)
                    
                    
                    
                    merged_msa_to_info [data['merged_msa']] = [pid, date, gene, compartment, replicate]
                    current_step = "Checking for merged_json"
                                    
                    files_to_copy = [data['merged_json']]
                    # read coverage info
                    
                    ci = None
                    
                    current_step = "Opening merged_json"
                    with open (files_to_copy[0]) as fh: 
                        current_step = "Loading %s" % files_to_copy[0]
                        ci = coverage_info(json.load (fh))
                        
                                                  
                    if ci[0] is not None:    
                        
                        if 'overall_region' in data:
                            files_to_copy.append (data['overall_region'][1])
                            spanned_region = "-".join([str (k) for k in ci[0]])
                            
                            cell_entry = {'text' : spanned_region,
                                          'link' : join(dir_name,split(data['overall_region'][1])[1]),
                                          'target' : dir_name + "_" + spanned_region + ".fas"}
                                          
                            #print (cell_entry)
                            
                            row.append (cell_entry)
                        else:
                            row.append ("-".join([str (k) for k in ci[0]]))
                           
                        row.append (ci[1]['median'])
                    
                        files_to_copy.append (data['tn93_json'])
                        with open (files_to_copy[-1]) as fh:
                            current_step = "Loading %s" % files_to_copy[-1]
                            tn93 = json.load (fh)
                        
                        tn93_dist = tn93['Mean distance']
                        row.append (tn93_dist*100.)    
                        
                        files_to_copy.append (data['region_merged_processed'])
                        
                        
                        with open (files_to_copy[-1]) as fh:
                            current_step = "Loading %s " %  files_to_copy[-1]
                            diversity = json.load (fh)
                                 
                    
                        if len (diversity) == 0:
                            row.append (None)
                            row.append (None)
                        else:       
                            row.append (diversity_info (diversity, 'S')['max']*100.)
                            row.append (diversity_info (diversity, 'NS')['max']*100.)
                        
                        if 'tropism' in data and data['tropism']:
                            current_step = "Copying tropism data"
                            row.append (data['tropism']['X4'] * 100.)
                            has_tropism = True
                        else:
                            row.append ('N/A')

                        files_to_copy.append (data['json_rates'])
                        with open (files_to_copy[-1]) as fh:
                            current_step = "Loading json_rates %s" % files_to_copy[-1]
                            json.load (fh)
                                 
                        #print (row)         
                        if intrahost_info is not None:
                            try:
                                this_record = intrahost_info[pid][gene]
                                if has_compartment:
                                    this_record = this_record [compartment]
        
                                if has_replicate:
                                    this_record = this_record[replicate]

                                
                                with open (this_record['rates.tsv'], 'r') as fh:
                                    intra_host = csv.reader (fh, delimiter = '\t')
                                    headers = next (intra_host)
                                    if len (headers) == 14:
                                        add_key_if_missing (processed ['intrahost'], pid, {})
                                        store_here = add_key_if_missing (processed ['intrahost'][pid], gene, {})
                            
                                        if has_compartment:
                                            store_here = add_key_if_missing (store_here, compartment, {})
                                        if has_replicate and not ignore_replicate:
                                            store_here = add_key_if_missing (store_here, replicate, {})
                        
                                        for r in intra_host:
                                            div_date = parse_a_date (r[0])
                                            store_here[div_date] = {}
                                            for c in range(1,len(r)):
                                                 store_here[div_date][headers[c]] = r[c]
                                    
                                        
                                        
                            except (KeyError, AttributeError, TypeError, FileNotFoundError, ValueError) as e:
                                #print ("Failing key %s (%s)" % (gene, str (e)), file = sys.stderr)
                                pass
                    else:
                        row.extend ([None,ci[1]['median'],None,None,None,None])
                    
                    
                except (KeyError, AttributeError, TypeError, ValueError) as e:
                    if current_step:
                        row.extend ([None for k in range (len (row), len (processed['columns'])-1)])
                        print ("Failing %s" % current_step, row, e)
                    else:
                        continue

                result_path = os.path.join (store_dir,dir_name)
                if not os.path.exists (result_path):
                    os.makedirs (result_path)
            
                for file in files_to_copy:
                    shutil.copy (file, result_path)
                    
                row.append (dir_name)
                #print (dir_name, row)
                processed['data'].append (row)
                    
        except (KeyError,TypeError) as e:
            if record != 'F_ST':
                print ('Missing required record fields for %s' % (record), e, file = sys.stderr)
           
    if not has_tropism:
        for row in processed['data']:
            row.pop (-2)
        processed['columns'].pop (-2)

    print ("Compiled data on %d out of %d present individual NGS libraries" % (len (processed['data']), tried), file = sys.stderr)
           
    if intrahost_info is not None and pairwise_distance_info is not None:
        for d in processed['data']:

            try:
                store_dict = processed ['intrahost'][d[0]][d[2]]
                id = 2
                if has_compartment:
                    id += 1
                    store_dict = store_dict[d[id]]
                if has_replicate and not ignore_replicate:
                    id += 1
                    store_dict = store_dict[d[id]]
                
                
                tn93 = d[id + 3]
                if tn93 is None:
                    continue

                store_dict[d[1]]['tn93_diversity'] = tn93 * 0.01
                #print (store_dict, "\n", d[1], "\n\n")
            except (KeyError, AttributeError, TypeError) as e:
                continue
            
        for pair, distance_info in pairwise_distance_info.items():
            if short:
                pair_info = pair.split ('-')
                if len (pair_info) == 3:
                    tag1 = id_to_pid_date[int (pair_info[0])]
                    tag2 = id_to_pid_date[int (pair_info[1])]
                    gene = pair_info[2]
                    offset = 0
                else:
                    continue
            else:
                pair_info = pair.split ('|')
                if len (pair_info) == 2:
                    tag1 = merged_msa_to_info [pair_info[0]]
                    tag2 = merged_msa_to_info [pair_info[1]]
                    gene = tag1 [2]
                    offset = 1
                else:
                    continue
                    
            
            pid = tag1[0]
                                                  
            #print (tag1, tag2)
                                        
            try:
                store_dict = processed ['intrahost'][pid][gene]
                if has_compartment:
                    #store_dict = store_dict[tag1[3]]
                    if tag1[3] != tag2[3]:
                        continue
                '''
                if has_replicate:
                    store_dict = store_dict[tag1[4]]
                    if tag1[4] != tag2[4]:
                        continue
                '''
            
                #print (tag1, tag2)
                #store_dict [earilest_date_by_pid[pid]]['tn93_divergence'] = 0.0
                
                if pid == tag2 [0]:
                    baseline_date = earilest_date_by_pid[pid]
                    if baseline_date == tag1[1] or baseline_date[pid] == tag2[1]:
                 
                        #add_key_if_missing (store_dict, tag1[1], {})
                        #add_key_if_missing (store_dict, tag2[1], {})
                        #store_dict [tag1[1]]['tn93_divergence'] = 0.
                        #store_dict [tag2[1]]['tn93_divergence'] = 0.
                        
                        if earilest_date_by_pid[pid] == tag1[1]:
                            store_tag = tag2
                        else:
                            store_tag = tag1
                            
                        store_here = store_tag [1]
                        if has_compartment:
                            store_dict = store_dict[store_tag[3]]
                        if has_replicate and not ignore_replicate:
                            store_dict = store_dict[store_tag[4]]
                        
                        
                        add_key_if_missing (store_dict, store_here, {})
                        store_dict [store_here]['tn93_divergence'] = distance_info['Mean']
                        
                        if 'Histogram' in distance_info:
                           store_dict [store_here]['tn93_divergence_histogram'] =  distance_info['Histogram']
                           
                        #26919/20010921/BP/2/env/   
                        '''if '26919' in processed ['intrahost']:
                            if 'env' in processed ['intrahost']['26919']:
                                if 'BP' in processed ['intrahost']['26919']['env']:
                                    if '2' in processed ['intrahost']['26919']['env']['BP']:
                                        if '2001/09/21' in processed ['intrahost']['26919']['env']['BP']['2']:
                                            print ('BARF', pair, tag1, tag2)
                                            sys.exit (1)
                        '''
                        
            except (KeyError, AttributeError, TypeError) as e:
                #print (pid, e)
                pass
                    
    if has_fst:
        for f_pid, by_date in cache_file["F_ST"].items():
            for f_date, by_gene in by_date.items():
                for f_gene, by_pair in by_gene.items():
                    for pair_key, pair_data in by_pair.items():
                        try:
                            id1, id2 = pair_key.split ('|')
                            info1 = merged_msa_to_info[id1]
                            info2 = merged_msa_to_info[id2]
                            store_here = add_key_if_missing (processed["F_ST"], f_pid, {})
                            store_here = add_key_if_missing (store_here, parse_a_date(f_date), {})
                            store_here = add_key_if_missing (store_here, f_gene, [])
                            store_here.append ([info1[-1], info2[-1], pair_data])
                            

                        except:
                            raise
                    
         
    processed['data'].sort (key = lambda row : row[0: (3 + (1 if has_compartment else 0) + (1 if has_replicate else 0))])
    json.dump (processed, out, indent=1)

                                   
    
    return 0

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='scan the directory of NGS result files'
    )
 
    parser.add_argument(
        '-o', '--output',
        metavar='DIR',
        type=str,
        help='the directory to store web-ready files to',
        required = True,
    ) 
    
    parser.add_argument (
        '-j', '--json',
        type=argparse.FileType('w'),
        help='write the JSON file suitable for plotting here',
        required = True,
    ) 
 
    parser.add_argument (
        '-c', '--cache',
        type=argparse.FileType('r'),
        help='the pipeline cache file',
        required = True,
    ) 

    parser.add_argument (
        '-d', '--diversity',
        type=argparse.FileType('r'),
        help='the JSON file with intra-host evolution information',
        required = False,
    ) 

    parser.add_argument (
        '-w', '--pairwise',
        type=argparse.FileType('r'),
        help='the JSON file with pairwise distance information',
        required = False,
    ) 

    parser.add_argument(
        '-p', '--compartment',
        action='store_true',
        help='does the directory structure include compartment information',
    )    

    parser.add_argument(
        '-r', '--replicate',
        action='store_true',
        help='does the directory structure include replicate information',
    )    

    parser.add_argument(
        '-s', '--short',
        help='assume that pairwise results are keyed on id-id-gene (otherwise path|path)',
        action = 'store_true',
        default=False
    )


    parser.add_argument(
        '-g', '--ignore_replicate',
        help='ignore replicate information when computing intra-host tables',
        action = 'store_true',
        default=False
    )

    args = None
    retcode = -1
    args = parser.parse_args()
        
    #loadDRM ("Scores_PI.txt", "PR", "PI")
    #loadDRM ("Scores_NRTI.txt", "RT", "NRTI")
    #loadDRM ("Scores_NNRTI.txt", "RT", "NNRTI")
        
    retcode = main(json.load (args.cache), args.json, args.output, args.compartment, args.replicate, json.load(args.diversity) if args.diversity is not None else None , json.load(args.pairwise) if args.pairwise is not None else None, args.short, args.ignore_replicate)
    
    sys.exit(retcode)


