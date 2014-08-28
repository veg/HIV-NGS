import os, argparse, re, sys
import subprocess, time, csv
import json, operator
import itertools, datetime

    
def describe_vector (vector):
    vector.sort()
    l = len (vector)
    if l:
        return {'count': l, 'min': vector[0], 'max': vector[-1], 'mean': sum(vector)/l, 'median':  vector [l//2] if l % 2 == 1 else 0.5*(vector[l//2-1]+vector[l//2]), "IQR": [vector [l//4], vector [(3*l)//4]] }
    else:
        return {'count': l, 'min': None, 'max': None, 'mean': None, 'median':  None, "IQR": [None,None] }
       
    
def main (ngs_cache, pairwise_cache, csv_writer, short_format):
    
    print ("Loaded cache info on %d NGS runs" % len (ngs_cache), file = sys.stderr)                    
    print ("Loaded cache info on %d pairwise comparisons" % len (pairwise_cache), file = sys.stderr)            
    
    id_to_path = {}
    ids_seen = set ()
    if short_format:
    
    
        for path,d in ngs_cache.items():
            if type (d) == dict:   
                id = int(d['id'])
                if id in ids_seen:
                    print ("BARF %d" % id)
                    sys.exit (0)
                ids_seen.add (id) 
                id_to_path [id] = path
    else:
        for path,d in ngs_cache.items():
            if type (d) == dict:   
                for k,v in d.items():
                    if type (v) == dict:
                        if 'merged_msa' in v:
                            id_to_path [v['merged_msa']] = (path, k)

    processed = {}
    
    #test_pid = set (['050102333', '050108085'])
    #test_runs = set ()
    
    tbyg = {'gag' : "1.5%", 'rt' : "1.5%", 'env' : "5%"}
    
    unique_pids = set ()
    threshold   = 0.06

    for pair,value in pairwise_cache.items():
        if value is not None and value["Pair Count"] > 10000:
            
            if short_format:
                id0, id1, gene = pair.split ('-')
                pids  = [ngs_cache[id_to_path[k]]['patient_id'] for k in [int(id0), int(id1)]]
                dates = [ngs_cache[id_to_path[k]]['sample_date'] for k in [int(id0), int(id1)]]
            else:
                paths = pair.split ('|')
                pids  = [ngs_cache[id_to_path[k][0]]['patient_id'] for k in paths]
                dates = [ngs_cache[id_to_path[k][0]]['sample_date'] for k in paths]
                gene = id_to_path[paths[0]][1]
            
            if pids[0] < pids[1]:
                tag = (pids[0], pids[1], dates[0], dates[1])
            elif pids [0] > pids[1]:
                tag = (pids[1], pids[0], dates[1], dates[0])
            else:
                tag = (pids[0],pids[0], min (dates), max(dates))
        
            unique_pids.update ((pids[0], pids[1]))
         
            if tag not in processed:
                processed[tag] = {}
                
            if gene not in  processed[tag]:
                processed[tag][gene] = []
        
            link_ratio = value[tbyg[gene]] / value ["Pair Count"]
            processed[tag][gene].append (link_ratio)
            #if link_ratio > 0.06:
            #    print (",".join(paths), link_ratio)
    

    replicates = {'rt' : [], 'gag': [], 'env' : []}
    key_gene = 'rt'
    all_calls   = {'rt' : [], 'gag': [], 'env' : []}
    
    for t,v in processed.items():
        if t[0] != t[1] and key_gene in v and max (v[key_gene]) >= threshold:
            ratio = {key_gene : len ([k for k in v[key_gene] if k >= threshold]) / len (v[key_gene]) }
            repl = len (v[key_gene]) > 1
            if repl:
                replicates[key_gene].append (ratio[key_gene])
            all_calls[key_gene].append (ratio[key_gene])
            for sg in ['env','gag','rt']:
                if sg == key_gene: continue   
                if sg in v:
                    ratio[sg] = (len ([k for k in v[sg] if k >= threshold]) / len (v[sg]))
                    if repl:
                        replicates[sg].append (ratio[sg])
                    all_calls[sg].append (ratio[sg])
                else:
                    ratio[sg] = None
                    all_calls[sg].append (None)
                        
            '''
            if csv_writer is not None and True in [True if (ratio[gene] is not None and ratio[gene] >= 0.5) else False for gene in ['rt','env','gag']]:
                csv_writer.writerow (['|'.join ([t[0],aeh_time (t[2])]), '|'.join ([t[1],aeh_time (t[3])]), 0.01499999])      
            '''   
            
            non_null = [k for k in ratio.values() if k is not None]
            if csv_writer is not None and len (non_null) > 2 and len ([k for k in non_null if k >= 0.5])>= 2:
                csv_writer.writerow (['|'.join ([t[0],aeh_time (t[2])]), '|'.join ([t[1],aeh_time (t[3])]), 0.01499999])      
            
            #print (ratio)
                
    
    print ("%d unique PIDs" % len (unique_pids))
    
    confirmed = {'or' : [0,0]}
    
    for i in range (len(all_calls[key_gene])):
        or_tag = 0
        for sg in ['env','gag','rt']:
            if sg == key_gene: continue   
            if sg not in confirmed:
                confirmed[sg]=[0,0]
                
            if all_calls [sg][i] is not None:
                confirmed[sg][1 if all_calls [sg][i] > 0 else 0] += 1
                if all_calls [sg][i] > 0:
                    or_tag = 1
            
        confirmed['or'][or_tag] += 1
        
    print (confirmed)
    print (describe_vector(replicates[key_gene]))
    
                
    return 0
    
def aeh_time (t):
    return datetime.datetime.strptime (t,"%Y%m%d").strftime ("%m%d%Y")

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='examine pairwise distance comparisons on NGS reads'
    )
    parser.add_argument(
        '-i', '--input',
        metavar='NGS',
        type=argparse.FileType ('r'),
        help='the .json cache of the NGS pipeline',
        required = True,
    )
    
    parser.add_argument(
        '-c', '--cache',
        metavar='JSON',
        type=argparse.FileType ('r'),
        help='the .json cache for the pairwise comparison script',
        required = True,
    )

    parser.add_argument(
        '-d', '--distances',
        metavar='CSV',
        type=argparse.FileType ('w'),
        help='the CSV cache for the pairwise distances',
        required = False,
    )
    
    parser.add_argument(
        '-s', '--short',
        help='store results in a dict keyed on id-id-gene (otherwise path-path)',
        action = 'store_true',
        default=False
    )
    
    args = None
    retcode = -1
    args = parser.parse_args()
    
    csv_writer = None
    if args.distances:
        csv_writer = csv.writer (args.distances)
        csv_writer.writerow (['ID1','ID2','Distance'])
        
    retcode = main(json.load (args.input), json.load (args.cache),csv_writer, args.short)
    
    sys.exit(retcode)


