import os, argparse, json, sys



if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='view (and optionally delete) a path key from a .json cache'
    )
 
    parser.add_argument (
        '-c', '--cache',
        type=argparse.FileType('r+'),
        help='the pipeline cache file',
        required = True,
    ) 

    parser.add_argument (
        '-p', '--path',
        type=str,
        help='the (partial) file path to retrieve',
        required = True,
    ) 

    parser.add_argument (
        '-g', '--gene',
        type=str,
        help='the gene to retrieve from the path dict',
        required = True,
    ) 


    parser.add_argument (
        '-d', '--delete',
        action='store_true',
        help='whether or not to delete the key',
    ) 

    args = None
    retcode = -1
    args = parser.parse_args()
        
    main_json = json.load (args.cache)
    
    matched_keys = set ()
    for key in main_json:
        if key.find (args.path) >= 0:
            matched_keys.add (key)
            
    if len (matched_keys) == 1:
        the_key = list (matched_keys)[0]
        if args.gene in main_json[the_key]:
            print (json.dumps (main_json[the_key][args.gene], indent = 2, sort_keys = True))
            if args.delete:
                args.cache.seek (0,0)
                del main_json[the_key][args.gene]
                args.cache.truncate ()        
                json.dump (main_json, args.cache, sort_keys=True, indent=2)
        else:
            print ("Gene key not found\n", json.dumps (main_json[the_key], indent = 2, sort_keys = True))
    else: 
        print ("Found multiple matching keys: ", matched_keys)
    
    retcode = 0
    
    sys.exit(retcode)


