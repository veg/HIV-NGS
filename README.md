HIV-NGS
=======

HIV NGS Processing pipeline

# 1 - Setup

## 1.1 - File Formats

Data prepared for this pipeline should be in paired .fna and .qual format. File pairs compressed into .zip format are acceptable.

## 1.2 -  Directory Structure

This pipeline is sensitive to the directory structure containing the data. Suppose all data files are contained in this directory:

    /data/NGS-DATA/

If the data files are not divided into compartments or replicates then each file pair should be contained in a sub-directory formatted as in this example:

    $ ls /data/NGS-DATA/PATIENT-ID/YYYYMMDD/COMPARTMENT/REPLICATE/
    example.fna example.qual

If compartment or replicate information about the data is not available, the appropriate level of sub-directory can be omitted.

The data directory should also contain a subfolder titled 'results', which should in turn contain a subfolder titled 'web-ready'. The following commands will create them:

    $ mkdir /data/NGS-DATA/results
    $ mkdir /data/NGS-Data/results/web-ready

# 2 - Processing

## 2.1 - Directory Scanner

The first step is to run directoryscanner.py on the directory containing the data. If the directory structure includes both compartment and replicate information, use the following command:

    $ cd /data/NGS-DATA/
    $ /opt/python-3.3.1/bin/python3 /opt/NGSpipeline/python/directoryscanner.py -i `pwd` -r `pwd` -c `pwd`/results/ds-cache.json -d ID DATE COMPARTMENT REPLICATE -n 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15

If the directory structure does not include replicate information, replace the value REPLICATE after the -d flag with NONE. If it does not include compartment information, replace the value COMPARTMENT after the -d flag with NONE.

## 2.2 - Intra Host Evolution

Once the data have been processed with directoryscanner.py, we run intra_host_evolution.py, still from the NGS-DATA directory

    $ /opt/python-3.3.1/bin/python3 /opt/NGSpipeline/python/intra_host_evolution.py -i `pwd`/results/ds-cache.json -c `pwd`/results/ihe-cache.json -p -d

Omit the -d flag if the directory structure does not include replicates and omit the -p flag if it does not include compartments.

## 2.3 - Pairwise Distances

Next we run pairwise_distances.py.

    $ /opt/python-3.3.1/bin/python3 /opt/NGSpipeline/python/pairwise_distances.py -i `pwd`/results/ds-cache.json -c `pwd`/results/pd-cache.json

# 3 - Postprocessing

Finally, we can use result_processor.py to produce web-ready output files.

    $ /opt/python-3.3.1/bin/python3 /opt/NGSpipeline/python/result_processor.py -o `pwd`/results/web-ready/ -j `pwd`/results/web-ready/results.json -c `pwd`/results/ds-cache.json -d `pwd`/results/ihe-cache.json -w `pwd`/results/pd-cache.json -p -r

Omit the -r flag if the directory structure does not include replicates and omit the -p flag if it does not include compartments.
