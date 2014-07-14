#!/usr/bin/env python3.2

import csv, argparse, sys, datetime, time, random, os.path, shutil, json

#-------------------------------------------------------------------------------		
    

arguments = argparse.ArgumentParser(description='Read filenames.')

arguments.add_argument('-i', '--input', help = 'Read UDS information from', required = True, type = argparse.FileType('r'))
arguments.add_argument('-l', '--lookup', help = 'Masterbook to AIEDRP PID', required = True, type = argparse.FileType('r'))
arguments.add_argument('-o', '--output', help = 'Send the output to ', required = True, type = str)
arguments.add_argument('-c', '--compartment', help = 'Default sampling compartment ', type = str, default = "BP")
arguments.add_argument('-j', '--json', help = '.json cache file', required = True, type = str)

settings = arguments.parse_args()

if settings.output == None:
	settings.output = sys.stdout


mb_to_aiedrp = {}

lookupReader = csv.reader(settings.lookup)
next(lookupReader)
for line in lookupReader:
    if len (line[2]):
        mb_to_aiedrp [line[2] if len (line[2]) == 4 else line[2][1:]] = line[1].replace('-','')
    else:
        mb_to_aiedrp [line[0]] = line[0]
        mb_to_aiedrp [line[0].upper()] = line[0]
  

udsReader = csv.reader(settings.input, delimiter = "\t")
    
paths_done = {}
base_path = "/data/sci/454data"

if os.path.exists(settings.json):
    with open (settings.json, "r") as fh:
        paths_done = json.load (fh)
else:
    paths_done = {}
        
if not os.path.exists (settings.output):
    os.makedirs (settings.output) 
    
for line in udsReader:
    for col_set in [[0,1,2], [3,4,5]]:
        file_name = os.path.normpath(os.path.join(base_path,line[col_set[0]]))
        if file_name not in paths_done:
            try:
                #check for the presence of [.fasta and .qual or a .fastq file]
                
                files_to_move = None
                for req_ext in [['.fna','.qual'],['.fastq']]:
                    file_exists = [os.path.exists(file_name+k) for k in req_ext]
                    if False not in file_exists:
                        files_to_move = [file_name+k for k in req_ext]
                        break
                    
                if files_to_move is not None:
                    pat_id = mb_to_aiedrp[line[col_set[1]].upper()]
                    sample_date = datetime.datetime.strptime(line[col_set[2]],"%Y-%m-%dT%H:%M:%S" ).strftime ("%Y%m%d")
                    replicate   = 1
                    result_path = os.path.join (settings.output, pat_id, sample_date, settings.compartment, str (replicate))
                    while os.path.exists (result_path):
                        replicate += 1
                        result_path = os.path.join (settings.output, pat_id, sample_date, settings.compartment, str (replicate))
                         
                    os.makedirs (result_path) 
                    
                    paths_done [file_name] = result_path
                    
                    for file in files_to_move:
                        name_to_move = os.path.split (file)
                        copy_to = os.path.join (result_path, name_to_move[1])
                        if not os.path.exists (copy_to):
                            os.symlink (file, copy_to)
                            print ("%s => %s" % (file, copy_to))
                    with open (settings.json, "w") as fh:
                        json.dump (paths_done, fh,sort_keys=True, indent=4)

            
                
            except KeyError as e:
                pass
            
