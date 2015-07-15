#!/usr/bin/env python3.2

import csv, argparse, sys, datetime, time, random, os.path, shutil, json,re

                    
##########
    
arguments = argparse.ArgumentParser(description='Make consensus from rates.json.')

def non_neg (string):
    r = float (string)
    if r < 0:
        raise argparse.ArgumentTypeError("'%s' is not a non-negative real" % string)
    return r

def valid_regexp (string):
    try:
        return re.compile (string)
    except:
        raise argparse.ArgumentTypeError("'%s' is not a valid regexp" % string)
        
    return None


def make_label (file, gene, tag):
    return '|'.join ([file['patient_id'], file['sample_date'], file['compartment'], file['replicate'], gene, tag])

def consensus (json_rates, cutoff, trim):
    sites = sorted([int (k) for k in json_rates['posteriors']])

    consensus  = []
    last_site = 0

    for s in sites:
        if s - last_site > 1:
            for k in range (last_site, s):
                consensus.append (['-',0])
        else:
            positions = [(k,v) for k,v in json_rates ['posteriors'][str(s)].items () if k != 'variants']
            max_c     = max (positions, key = lambda x : x[1])
            coverage  = sum ([k[1][2] for k in positions])
            consensus.append ([max_c[0], coverage])
        
        last_site = s
    
    if cutoff < 1:
        coverages = sorted([k[1] for k in consensus if k[1] > 0])
        cutoff = max (250, cutoff * (coverages[len (coverages) // 2] if len (coverages) % 2 else (coverages[len (coverages) // 2 -1] + coverages[len (coverages) // 2])/2. ))
    
    consensus = ''.join ([k[0] if k[1] >= cutoff else '-' for k in consensus])
    if trim:
        consensus = consensus.lstrip ("-").rstrip ("-")
        
    return consensus

if __name__ == '__main__':
  
    arguments.add_argument('-j', '--json', help = 'Main pipeline .json', required = True, type = argparse.FileType('r'))
    arguments.add_argument('-c', '--cutoff', help = 'punch positions with this coverage or worse: >1 absolute coverage; 0-1 : %% of median value', type = float, default = 0.25)
    arguments.add_argument('-t', '--truncate', help = 'strip leading and trailing gaps', action = "store_true", default = False)
    arguments.add_argument('-n', '--id', help = 'id filter', type = valid_regexp, default = None)
    arguments.add_argument('-g', '--gene', help = 'gene filter', type = valid_regexp, default = None)
    arguments.add_argument('-r', '--replicate', help = 'replicate filter', type = valid_regexp, default = None)
    arguments.add_argument('-p', '--compartment', help = 'comparment filter', type = valid_regexp, default = None)
    arguments.add_argument('-m', '--min_length', help = 'minimum number of non-gap chars', type = int, default = 50)
    arguments.add_argument('-o', '--clones', help = 'extract top X clones', type = int, default = 0)


    settings = arguments.parse_args()
  
  
    analysis_json = json.load (settings.json)

    reg_exp_filter = {'patient_id' : settings.id, 'replicate': settings.replicate, 'compartment': settings.compartment}

    for k,v in analysis_json.items():
        filter_pass = True
    
        for validator, filter in reg_exp_filter.items():
            if validator in v and filter and filter.search (v[validator]) is None:
                filter_pass = False
                break
            
        if filter_pass:
            for gene,data in v.items():
                try:
                    if type (data) is dict:
                            if settings.gene is None or settings.gene.search (gene):
                                if settings.clones > 0:
                                    if 'overall_region' in data and data['overall_region'] is not None:
                                        label = make_label (v, gene, "")
                                        try:
                                            has_merged_counts = data['merged_counts']['total']
                                        except:
                                            has_merged_counts = None
                                        with open (data['overall_region'][1], 'r') as fh:
                                            last_line = 0
                                            seq_data = []
                                            for line in fh:
                                                if last_line % 2 == 0:
                                                    tag = line.replace ("cluster_", label)
                                                    counts = line.split (':')
                                                    tag = counts[0].replace ("cluster_", label)
                                                    frac = 1
                                                    if len (counts) > 1:
                                                        frac = int (counts[1])
                                                        #tag += '|' + (counts[1] if has_merged_counts is None else "%.2f" % (100*float (counts[1]) / has_merged_counts))
                                                    seq_data.append ([tag,frac,None])
                                                else:
                                                    seq_data[-1][2] = line.rstrip()
                                                last_line += 1
                                        
                                            seq_data.sort (key = lambda a : -a[1])
                                            for seq in seq_data[0:settings.clones]:
                                                tag = seq[0] + '|' + (seq[1] if has_merged_counts is None else "%.2f" % (100*float (seq[1]) / has_merged_counts))
                                                print ("%s\n%s"%(tag,seq[2]))
                                else:
                                    if 'json_rates' in data and data['json_rates'] is not None:
                                        with open (data['json_rates'], 'r') as fh:
                                            c_seq = consensus (json.load (fh), settings.cutoff, settings.truncate)
                                            no_gap = c_seq.replace ('-','')
                                            if len (no_gap) > settings.min_length and len (c_seq) - len (no_gap) < settings.min_length*0.3:
                                                print (">%s\n%s" % (make_label (v, gene, "consensus"), c_seq))
                                
                except:
                    raise
                    pass
        
    
    
    