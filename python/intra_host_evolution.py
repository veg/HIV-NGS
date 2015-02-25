import os, argparse, re, sys
import subprocess, time, csv
import json, operator
import itertools
import queue
import threading
import hashlib
import ngs_task_runners as ntr

global task_queue
global host_cache

write_to_cache_lock = threading.Lock()
path_to_this_file = os.path.dirname(os.path.realpath(__file__))


def realign_the_files (node_to_run_on):
    global host_cache
    global task_queue

    while True:
        try:
            in_msa, reference, out_bam, out_msa, store_here, sample_date, file_hash, cache_path = task_queue.get()

            print ("Running bealign on %s (node %s)" % (in_msa, node_to_run_on), file = sys.stderr)
            
            with open (reference, 'r') as fh:
                lines = [k for k in fh.readlines()]
            
            with open (in_msa, 'r') as fh:
                lines.extend ([k for k in fh.readlines()])
                
            with open (out_msa, 'w') as fh:
                fh.writelines (lines)
                        
            command = ['/usr/bin/bpsh', str(node_to_run_on), '/opt/share/pytools/utils/python3.3/bealign', '-r', reference, '-e', '0.5', '-m', 'HIV_BETWEEN_F', out_msa, out_bam ]      
            subprocess.check_call (command, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL) 
            command = ['/opt/share/python3.3/bam2msa', out_bam, out_msa]
            subprocess.check_call (command,stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL) 


            write_to_cache_lock.acquire()
            store_here[in_msa]['hash'] = file_hash
            store_here[in_msa]['reference'] = reference
            store_here[in_msa]['bam'] = out_bam
            store_here[in_msa]['msa'] = out_msa
            store_here[in_msa]['sample_date'] = sample_date
            print ("Spooling updated JSON file with %d records" % len (host_cache), file = sys.stderr)
            refresh_json_file (cache_path, host_cache)
            write_to_cache_lock.release()
            
        except subprocess.CalledProcessError as e:
            print ('ERROR: process call failed',' '.join (command),file = sys.stderr)
               
        task_queue.task_done()

'''
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
'''


def run_hyphy_div (node_to_run_on):
    global host_cache
    global task_queue
    
    hbf = os.path.join (path_to_this_file, "../hyphy", "intra_host.bf")


    while True:
        try:
            files, cache_path = task_queue.get()
                
            #print ([[k, type (k)] for k in files.values() if type (k) == dict])
            alignment_files = sorted ([[k['msa'],k['sample_date']] for k in files.values() if type (k) == dict], key = lambda x : x[1])
            out_file = os.path.join (os.path.dirname (alignment_files[0][0]), 'rates.tsv')
            
            infile = "\n".join ([i for m in alignment_files for i in m] + [out_file, "EOF"])
            
            #print (infile)
            
            print ("Running diversity estimates for %s (node %s)" % (out_file, node_to_run_on), file = sys.stderr)
            command =  ['/usr/bin/bpsh', str(node_to_run_on), '/usr/local/bin/HYPHYMP', hbf]
            process = subprocess.Popen (command, stdin = subprocess.PIPE, stderr = subprocess.PIPE, stdout = subprocess.PIPE) 
            out, err = process.communicate (bytes(infile,'UTF-8'))
            err = err.decode('UTF-8')
            if len (err) and err.find ("Error") >= 0:
                print ("Call to %s with input '%s' generated an error: `%s`" % (' '.join (command), infile, err))
            
            write_to_cache_lock.acquire()
            files['rates.tsv'] = out_file
            print ("Spooling updated JSON file with %d records" % len (host_cache), file = sys.stderr)
            refresh_json_file (cache_path, host_cache)
            write_to_cache_lock.release()
             
        except subprocess.CalledProcessError as e:
            print ('ERROR: process call failed',' '.join (command),file = sys.stderr)
               
        task_queue.task_done()
                
def check_missing_key_file (values, key):
    try:
        return not os.path.exists(values[key])
    except:
        pass
    return True
        



def refresh_json_file (path, object):
    with open (path, "w") as fh:
        json.dump (object, fh, sort_keys=True, indent=1)
    
def drill_down (d, redo_overall, task_queue, host_cache_path):
    for k, v in d.items():
        if type (v) == dict:
            for k2 in v.values():
                if type (k2) == dict:
                    for k3,v3 in k2.items():
                        if k3 == 'rates.tsv':
                            continue
                        if type (v3) == dict:
                            drill_down (v, redo_overall, task_queue, host_cache_path)
                        elif type (v3) == str:
                            if set(v) in redo_overall or check_missing_key_file (v, 'rates.tsv'):
                                #print (v)
                                task_queue.put ([v, host_cache_path])
                            else:
                                redo_me = False
                                with open (v['rates.tsv']) as fh:
                                    lines = fh.readlines()
                                    if len (lines) != len (v):
                                        redo_me = True
                                        #print (v, lines)
                                if redo_me:
                                    print ('Redoing %s because it has the wrong number of lines' % v['rates.tsv'], file = sys.stderr)
                                    task_queue.put ([v, host_cache_path])    
                        break
                break
                    
      
    
def main (ngs_cache, host_cache_path, nodes_to_run_on, genes, has_compartment, has_replicate, file_key):
    
    print ("Loaded cache info on %d NGS runs" % len (ngs_cache), file = sys.stderr)            
             
    global host_cache
    global task_queue
    
    if os.path.exists(host_cache_path):
        with open (host_cache_path, "r") as fh:
            host_cache = json.load (fh)
    else:
        host_cache = {}
        
    print ("Loaded cache info on %d hosts" % len (host_cache), file = sys.stderr) 
    
   
    task_queue = queue.Queue ()
    for node in nodes_to_run_on:
        t = threading.Thread (target = realign_the_files, args = (node,))
        t.daemon = True
        t.start ()
     
        
    redo_overall = set()
        
    for gene in genes:
        print ("Scheduling %s file conversions..." % gene, file = sys.stderr)
        
        for id, status in ngs_cache.items():
            try:
                pid             = status['patient_id']
                sample_date     = status['sample_date']
                overall_path    = status[gene][file_key][1]
                reference_file  = status[gene]['reference_sequence']
                if has_compartment:
                    compartment = status['compartment']
                if has_replicate:
                    replicate = status ['replicate']
                    
                redo_this_file = False

                with open (overall_path, 'rb') as fh:
                    m = hashlib.md5()
                    m.update (fh.read())
                    file_hash = m.hexdigest()
 
                if pid not in host_cache:
                    host_cache[pid] = {}
                
                if gene not in host_cache[pid]:
                    host_cache[pid][gene] = {}
                    
                this_file_record =  host_cache[pid][gene]
                
                if has_compartment:
                    if compartment not in this_file_record:
                        this_file_record [compartment] = {}
                    this_file_record = this_file_record[compartment]
                
                if has_replicate:
                    if replicate not in this_file_record:
                        this_file_record[replicate] = {}
                    this_file_record = this_file_record[replicate]
                                        
                if overall_path not in this_file_record:
                    this_file_record[overall_path] = {}
                    
                    
                if check_missing_key_file (this_file_record[overall_path], 'bam') or check_missing_key_file (this_file_record[overall_path], 'msa') or file_hash != this_file_record[overall_path]['hash'] or sample_date != this_file_record[overall_path]['sample_date'] or reference_file != this_file_record[overall_path]['reference']: 
                    redo_this_file = True
                
                if redo_this_file:
                    task_queue.put ([overall_path, reference_file, overall_path + ".bam", overall_path + ".msa", this_file_record, sample_date, file_hash, host_cache_path])
                    redo_overall.update (set(this_file_record))
                    
            except (KeyError,TypeError) as e:
                pass    

            except Exception as e:
                print (id)
                raise (e)
                
    task_queue.join()   
    #print (redo_overall)
    
    refresh_json_file (host_cache_path, host_cache)
    
    task_queue = queue.Queue ()
    for node in nodes_to_run_on:
        t = threading.Thread (target = run_hyphy_div, args = (node,))
        t.daemon = True
        t.start ()

    drill_down (host_cache, redo_overall, task_queue, host_cache_path)
    
    
    task_queue.join()    
    refresh_json_file (host_cache_path, host_cache)

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
        help='the .json cache for the intrahost evolution script',
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
        '-n', '--node',
        type=str,
        help='run MP scripts on this node',
        default = "2"
    )        
    args = None
    retcode = -1
    args = parser.parse_args()
    
    retcode = main(json.load (args.input), args.cache, args.node.split (","), args.genes.split (','), args.compartment, args.replicate,"overall_region")
    
    sys.exit(retcode)


