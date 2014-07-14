/* check compartmentalization */

_in_directory = "/home/sergei/Projects/RonEllis454/results/NM487/20020128/";
genes           = {{"env","gag","rt"}};
suffixes        = {{"BP","CSF"}};
_scratch_file   = PATH_TO_CURRENT_BF + "_f_st_workfile";

_size   = 210;
_stride = 30;
_max    = 1000;

for (_g = 0; _g < Columns (genes); _g += 1) {
    _from = 0;
    _to = _from + _size;
    
    while (_to < _max) {
        
        for (_i = 0; _i < 2; _i += 1) {
            _in_file = _in_directory + DIRECTORY_SEPARATOR + suffixes[_i] + 
                        DIRECTORY_SEPARATOR + genes[_g] + DIRECTORY_SEPARATOR + "region_reduced_" + _from + "-" + _to + ".msa";   
            if ( (!_in_file) == 0) {
                break;
            }
            
            if (_i == 0) {
                DataSet seqs = ReadDataFile (_in_file);
            } else {
                DataSet second = ReadDataFile (_in_file);
            }
        }
        
        if (_i == 2) {
            if (seqs.species <= 3 || second.species <= 3) {
                break;
            }
            DataSet joint = Combine (seqs, second);
            for (_s = 0; _s < joint.species; _s += 1) {
                GetString (seqName, joint, _s);
                if (_s < seqs.species) {
                    seqName = suffixes[0] + "_" + seqName;
                } else {
                    seqName = suffixes[1] + "_" + seqName;
                }
                SetParameter (joint, _s, seqName);
            }
            DataSetFilter jointF = CreateFilter (joint,1);
            fprintf (stdout, genes[_g], ":", _from, "-", _to, ":", seqs.species, "/", second.species, "\n");
            fprintf (_scratch_file, CLEAR_FILE, jointF);
                            
            GLOBAL_FPRINTF_REDIRECT = "/dev/null";            
            
            ExecuteAFile (HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "F_ST.bf",
                          {"00" : "Distance formulae",
                         "01" : "Nucleotide/Protein",
                         "02" : _scratch_file,
                         "03" : "^" + suffixes[0],
                         "04" : "^" + suffixes[1],
                         "05" : "y",
                         "06" : "TN93",
                         "07" : "Skip",
                         "08" : "But of course",
                         "09" : "100"});
                 
            GLOBAL_FPRINTF_REDIRECT = None;            
            fprintf (stdout, resultAVL["P (Hudson, Slatkin and Maddison)"], "\n");          
        }
        
        _from += _stride;
        _to += _stride;
    }
}