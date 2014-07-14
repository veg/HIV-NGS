import os, argparse, re, sys, datetime
import subprocess, time, csv
import datetime, shutil, json, operator
from os.path import join, splitext, isfile
import math

genes = {'env', 'rt', 'gag'}

gene_ranges = {'rt' : range (54,232), 'gag': range (61,144), 'env': range (40,179)}

def describe_vector (vector):
    vector.sort()
    l = len (vector)
    return {'count': l, 'min': vector[0], 'max': vector[-1], 'mean': sum(vector)/l, 'median':  vector [l//2] if l % 2 == 1 else 0.5*(vector[l//2-1]+vector[l//2]), "IQR": [vector [l//4], vector [(3*l)//4]] }

def main (json_in, edi_info, dual_infection_info):
    
    median_coverage_by_gene = {}
    for g in genes:
        median_coverage_by_gene[g] = []
        
    info_by_patient = {}
    
    '''
    
    [id] : [date: ['gene': [replicate coverage],...]], ...]
    
    '''
    
    for path, value in  json_in.items():
        if type (value) != dict : continue
        pid = value['patient_id']
        sample_date = datetime.datetime.strptime(value['sample_date'], "%Y%m%d")
        
        if edi_info is not None:
            if pid in edi_info and 'EDI' in edi_info[pid]:
                sample_date = sample_date - datetime.datetime.strptime(edi_info[pid]['EDI'], "%Y-%m-%d")
                #print (type (sample_date) == datetime.timedelta)
            
        
        if pid not in info_by_patient:
            info_by_patient[pid] = {}
        if sample_date not in info_by_patient[pid]:
            info_by_patient[pid][sample_date] = {}
            
        
        for gene in genes:
            if gene in value and 'merged_json' in value[gene]:
                with open (value[gene]['merged_json'],'r') as fh:
                    if gene not in info_by_patient[pid][sample_date]:
                        info_by_patient[pid][sample_date][gene] = []
                        
                    coverage_json = json.load (fh)
                    median_coverage = []
                    for k in gene_ranges [gene]:
                        if str(k) in coverage_json:
                            median_coverage.append (sum(coverage_json[str(k)].values()))
                        else:
                            median_coverage.append (0)
                            
                    median_coverage = describe_vector (median_coverage) ['median']
                    median_coverage_by_gene[gene].append (median_coverage)
                    if median_coverage > 500:
                        info_by_patient[pid][sample_date][gene].append (median_coverage)
                    
    for g in genes:
        desc = describe_vector (median_coverage_by_gene[g])
        print ("%s, median %g, IQR %g-%g" % (g, desc['median'], desc['IQR'][0], desc['IQR'][1]))
   
    pids_seen = set ()   
    print ("\t".join (["PID (DI)","Sample date", "rt", "gag", "env"]))
    to_print =  {}
    for pid in info_by_patient:
        pids_seen.add (pid)
        #pid_label = str (len (pids_seen)) + ((" (%s)" % dual_infection_info[pid]) if dual_infection_info is not None and pid in dual_infection_info else "")
        pid_label = (("(%s) " % dual_infection_info[pid]) if dual_infection_info is not None and pid in dual_infection_info else "") + str (pid)
        pinfo = sorted (info_by_patient [pid].keys())
        to_print [pid_label] = []
        for i, d in enumerate(pinfo):
            rec = info_by_patient [pid][d]
            to_print[pid_label].append ('\t'.join ([pid_label if i == 0 else '', str (d.days // 7) + " wks post EDI" if type (d) == datetime.timedelta else d.strftime ("%m/%Y")] +
                            [(str (len (rec[gene])) + " (" + str(describe_vector (rec[gene])['median']) + ")") if gene in rec and len (rec[gene]) > 0 else "0" for gene in ['rt','gag','env']
                        ]))
    
    for a_str in sorted (to_print.keys()):
        for line in to_print[a_str]:
            print (line)
    return 0

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='scan the directory of NGS result files'
    )
 
    parser.add_argument(
        '-i', '--input',
        metavar='DIR',
        type=argparse.FileType ('r'),
        help='the .json cache of the NGS pipeline',
        required = True,
    )
    
    parser.add_argument(
        '-e', '--edi',
        metavar='DIR',
        type=argparse.FileType ('r'),
        help='the .json file with EDI information',
        required = False,
    )    

    parser.add_argument(
        '-d', '--dual',
        metavar='DIR',
        type=argparse.FileType ('r'),
        help='the .csv file with DI information',
        required = False,
    )    
 
    args = None
    retcode = -1
    args = parser.parse_args()
    
    di = None
    if args.dual is not None:
        di_reader = csv.reader (args.dual)
        next (di_reader)
        di = {}
        for l in di_reader:
            di[l[0]] = l[1]
        
    retcode = main(json.load (args.input), json.load (args.edi) if args.edi is not None else None, di)
    
    sys.exit(retcode)


