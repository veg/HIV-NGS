import os, argparse, re, sys
import subprocess, time, csv
import json, operator
import itertools
import queue
import threading

global task_queue
global pairwise_cache

write_to_cache_lock = threading.Lock()
last_spool = time.time()    

def_genes = {'env' : 0.03, 'gag' : 0.015, 'rt' : 0.015, 'pr': 0.015}
  
def pairwise_distance_estimator (node_to_run_on):
    global last_spool

    while True:
        pairwise_key, in_path1, in_path2, gene, threshold, pairwise_cache_path, store_histogram,  subsample = task_queue.get()

        if threshold is None:
            threshold =  def_genes[gene] if gene in def_genes else 0.015;

        print ("Running tn93 on %s and %s [%s] with threshold = %g using node %s subsampling rate % g (queue size = %d)" % (in_path1, in_path2, gene, threshold, node_to_run_on, subsample, task_queue.qsize()), file = sys.stderr)
        try:
            process = subprocess.Popen (['/usr/bin/bpsh', node_to_run_on, '/usr/local/bin/tn93', '-t', str(threshold), '-l', '200', '-u', str(subsample), '-c', '-q', '-s', in_path2, in_path1], stdin = subprocess.DEVNULL, stderr = subprocess.PIPE, stdout = subprocess.PIPE, universal_newlines = True) 
            out, json_out = process.communicate ()
            full_json = json.loads(json_out)
            baseline_json = {}
            baseline_json ["Pair Count"] = full_json ["Comparisons accounting for copy numbers "] 
            baseline_json ["Mean"] = full_json ["Mean distance"]
            for upper_bound in range (1,11):
                baseline_json ["%g%%" % (upper_bound * 0.5)] = sum ([k[1] for k in full_json["Histogram"][:upper_bound]])
            if store_histogram:
                baseline_json ["Histogram"] = full_json["Histogram"]
            
            
        
        except subprocess.CalledProcessError as e:
            print ('ERROR: tn93 call failed',e,file = sys.stderr)
            baseline_json = None
        except ValueError as e:
            baseline_json = None
        
        write_to_cache_lock.acquire()
        pairwise_cache [pairwise_key] = baseline_json
        current_time = time.time()
        if current_time - last_spool > 300:
            print ("Spooling updated JSON file with %d records" % len (pairwise_cache), file = sys.stderr)
            refresh_json_file (pairwise_cache_path, pairwise_cache)
            last_spool = current_time
            
        write_to_cache_lock.release()
                
        task_queue.task_done()


def refresh_json_file (path, object):
    with open (path, "w") as fh:
        json.dump (object, fh, sort_keys=True, indent=2)
    
       
    
def main (ngs_cache, pairwise_cache_path, nodes_to_run_on, genes, file_key, intrahost_only, do_short, store_histogram, subsample):
    
    print ("Loaded cache info on %d NGS runs" % len (ngs_cache), file = sys.stderr)            
             
    global pairwise_cache
                
    if os.path.exists(pairwise_cache_path):
        with open (pairwise_cache_path, "r") as fh:
            pairwise_cache = json.load (fh)
    else:
        pairwise_cache = {}
        
    print ("Loaded cache info on %d pairwise comparisons" % len (pairwise_cache), file = sys.stderr)            
        
    
    for node in nodes_to_run_on:
        t = threading.Thread (target = pairwise_distance_estimator, args = (node,))
        t.daemon = True
        t.start ()
        
    for gene in genes:
        print ("Scheduling %s comparisons..." % gene, file = sys.stderr)
        
        for pair in itertools.combinations (ngs_cache.keys(), 2):
            try:
                file1 = ngs_cache[pair[0]][gene][file_key]               
                file2 = ngs_cache[pair[1]][gene][file_key]         
                
                id1 = ngs_cache[pair[0]]['id'] if do_short else file1
                id2 = ngs_cache[pair[1]]['id'] if do_short else file2 
                
                if intrahost_only:
                    try:
                        if ngs_cache[pair[0]]['patient_id'] !=  ngs_cache[pair[1]]['patient_id']:
                            continue
                    except Exception as e:
                        raise (e)
                
                if do_short:
                    if id1 < id2:
                        pairwise_key = "%d-%d-%s" % (id1, id2, gene)
                    else:
                        pairwise_key = "%d-%d-%s" % (id2, id1, gene)
                else:
                    if id1 < id2:
                        pairwise_key = "%s|%s" % (id1, id2)
                    else:
                        pairwise_key = "%s|%s" % (id2, id1)
                    
                    
                if pairwise_key not in pairwise_cache:
                    task_queue.put ([pairwise_key, file1, file2, gene, None, pairwise_cache_path, store_histogram, subsample])
                    #pairwise_cache [pairwise_key] = pairwise_distance_estimator (file1, file2, gene)
                            
                
            except (KeyError,TypeError):
                pass    
            except Exception as e:
                print (pair, gene, file_key)
                raise (e)
    
    task_queue.join()    
    refresh_json_file (pairwise_cache_path, pairwise_cache)
    return 0

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='run pairwise distance comparisons on NGS reads'
    )
    parser.add_argument(
        '-i', '--input',
        metavar='DIR',
        type=argparse.FileType ('r'),
        help='the .json cache of the NGS pipeline',
        required = True,
    )
    
    parser.add_argument(
        '-c', '--cache',
        metavar='JSON',
        type=str,
        help='the .json cache for the pairwise comparison script',
        required = True,
    )
    
    parser.add_argument(
        '-g', '--genes',
        metavar='genes',
        type=str,
        help='the comma separated list of genes to include in the comparison',
        default='env,gag,rt'
    )

    parser.add_argument(
        '-m', '--msa',
        metavar='MSA',
        type=str,
        help='the name of the merged multiple sequence alignment file to run the distance comparisons on',
        default='merged_msa'
    )
    
    parser.add_argument(
        '-t', '--intrahost',
        help='only perform distance calculations on intra-host samples',
        required = False, 
        action = 'store_true',
        default = False
    )
 
    parser.add_argument(
        '-s', '--short',
        help='store results in a dict keyed on id-id-gene (otherwise path-path)',
        action = 'store_true',
        default=False
    )
    
    parser.add_argument(
        '-u', '--subsample',
        type=float,
        help='subsample sequences with the probability',
        default = "1."
    )

    parser.add_argument(
        '-r', '--histogram',
        help='store the full histogram in the JSON',
        action = 'store_true',
        default=False
    )    
#    arguments.add_argument('-j', '--json', help = 'Output the network report as a JSON object', required = False,  action = 'store_true', default = False)
    parser.add_argument(
        '-n', '--node',
        type=str,
        help='run MP scripts on this node',
        default = "2"
    )        
    args = None
    retcode = -1
    args = parser.parse_args()
    
    task_queue = queue.Queue ()
    retcode = main(json.load (args.input), args.cache, args.node.split (","), args.genes.split (','), args.msa, args.intrahost, args.short, args.histogram, args.subsample)
    
    sys.exit(retcode)


