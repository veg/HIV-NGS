MESSAGE_LOGGING = 0;
LF_SMOOTHING_SCALER    = 0.1;
OPTIMIZATION_PRECISION = 0.01;

LoadFunctionLibrary ("HXB2Mapper", {"0": "Universal"});
LoadFunctionLibrary ("AncestralMapper");
LoadFunctionLibrary ("LocalMGREV");
LoadFunctionLibrary ("CF3x4");
LoadFunctionLibrary ("NJ");
LoadFunctionLibrary ("TreeTools");
LoadFunctionLibrary ("TreeFunctions");
LoadFunctionLibrary ("dSdNTreeTools");
LoadFunctionLibrary ("DescriptiveStatistics");
LoadFunctionLibrary ("BranchLengthFitters");
LoadFunctionLibrary ("CodonTools");
LoadFunctionLibrary ("GrabBag");
LoadFunctionLibrary ("_cwf.ibf");


file_list = {};
count_by_date = {};

while (1) {
    fscanf (stdin, "String", file);
    fscanf (stdin, "String", date_string);
    if (date_string == "EOF") {
        output_file = file;
        break;
    } else {
        ExecuteAFile (HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "CleanStopCodons.bf",
                {"1" : file,
                 "0" : "Universal",
                 "2" : "No/No",
                 "3" : file});    
    }
    
    count_by_date [date_string] += 1;
    
    if (Abs (file_list [date_string] ) == 0) {
        file_list [date_string] = file;
    } else {
        existing_file = file_list [date_string];
        DataSet existingData = ReadDataFile (existing_file);
        exclude = existingData.species;
        DataSet newData = ReadDataFile (file);
        
        
        for (k = 0; k < newData.species; k+=1) {
            GetString (sn, newData, k);
            SetParameter (newData, k, "set" +  count_by_date[date_string] + "_" + sn);
        }
        
        DataSet existingData = Combine (existingData, newData);
        
        
        DataSetFilter merged = CreateFilter (existingData,1,"",speciesIndex != exclude);
        if ((existing_file$"_merged$")[0] < 0) {
            existing_file = existing_file + "_merged";
        }
        fprintf (existing_file, CLEAR_FILE, merged);
        file_list [date_string] = existing_file;
        
    }
} 

unsortedDates = Rows (file_list);
uniqueDates = _sortStrings (unsortedDates);

//fprintf (stdout, file_list, "\n", uniqueDates, "\n");

excludeThese = {"0" : 1};

total_done = 0;

for (r = 0; r < Abs (file_list); r += 1) {
    if (r == 0) {
        DataSet allData = ReadDataFile (file_list[uniqueDates[r]]);
        for (s = 1; s < allData.species; s += 1) {
            GetString (sn, allData, s);
            SetParameter (allData, s, unsortedDates[r] + "_" + sn); 
        }
    } else {
        DataSet thisFile = ReadDataFile (file_list[uniqueDates[r]]);
        excludeThese [allData.species] = 1;
        for (s = 1; s < thisFile.species; s += 1) {
            GetString (sn, thisFile, s);
            SetParameter (thisFile, s, unsortedDates[r] + "_" + sn); 
        }
        DataSet allData = Combine (allData, thisFile);
    }
}

DataSetFilter allDataCodon = CreateFilter (allData,3,"",excludeThese[speciesIndex] == 0,GeneticCodeExclusions);
DataSetFilter filteredData = CreateFilter (allData,1,"",excludeThese[speciesIndex] == 0);

_normalizeSequenceNamesInAFilter ("filteredData");
allTree = InferTreeTopology (1);
Tree allTree = allTree;
allTreeAVL = allTree ^ 0;
bls = getNucRevBranchLengthsAndParameters ("filteredData", "allTree");

fprintf (output_file + ".nwk", CLEAR_FILE, PostOrderAVL2StringDistances (allTreeAVL, bls["lengths"]));

GetDataInfo (site_map, allDataCodon);

all_gap_map = -1;

for (s = 0; s < Columns (site_map); s += 1) {
    for (q = 0; q < allDataCodon.species; q += 1) {
        GetDataInfo (sq, allDataCodon, q, site_map[s]);
        if (+sq != Rows (sq)) {
            break;
        }
    }
    if (q == allDataCodon.species) {
        all_gap_map = site_map[s];
        break;
    }
}

global_partititon_string = "";
if (all_gap_map >= 0) {
    global_partititon_string * 128;
    for (s = 0; s < Columns (site_map); s += 1) {
        if (site_map[s] != all_gap_map) {
            global_partititon_string * ("," +  s*3);
            global_partititon_string * ("," +  (1+s*3));
            global_partititon_string * ("," +  (2+s*3));
        }
    }
    global_partititon_string * 0;
    global_partititon_string = global_partititon_string[1][Abs(global_partititon_string)-1];
} 

COUNT_GAPS_IN_FREQUENCIES = 0;

HarvestFrequencies (positionalFrequencies, allData,3,1,1);

_blStencils = ComputeScalingStencils(0);


nuc3x4 = CF3x4 (positionalFrequencies,GeneticCodeExclusions);
PopulateModelMatrix ("MGLocalQ", nuc3x4);
vectorOfFrequencies = BuildCodonFrequencies(nuc3x4);
Model MGLocal = (MGLocalQ, vectorOfFrequencies,0);

include_sequences = {};
for (k = 0; k < allData.species; k += 1) {
    if (excludeThese[k] == 0) {
        include_sequences + k;
    }
}

headers = {{"Date", "s_diversity", "ns_diversity",  "total_diversity",
                    "ds_diversity", "dn_diversity", 
                    "Length", "PNGS", "IsoelectricPoint",
                    "s_divergence","ns_divergence","total_divergence",
                    "ds_divergence", "dn_divergence"}};


call_count = 0;

for (date_index = 0; date_index < Rows (uniqueDates); date_index += 1) {
    thisDate = uniqueDates[date_index];
    
    DataSet thisFile = ReadDataFile (file_list[thisDate]);
        
    if (thisFile.species > 2) {
        if (date_index == 0) {
            reroot = None;
            GetString (reroot, thisFile, 0);
            rootSeq = makePartitionBuiltTreeFitMG ("thisFile", global_partititon_string, 
            speciesString,
            "sampled`thisDate`", reroot);
            fprintf (output_file, CLEAR_FILE, Join ("\t", headers), "\n");
        }    
    
        timepoint_info = makePartitionBuiltTreeFitMG ("thisFile", global_partititon_string, 
            "1-" + (thisFile.species),
            "sampled`thisDate`", None);
        
        call_count += 1;
        div_data = divergenceForFilter ("sampled`thisDate`_filter", rootSeq);
    } else {
        if (date_index == 0) {
            GetString (rootSeq, thisFile, 0);
            fprintf (output_file, CLEAR_FILE, Join ("\t", headers), "\n");
        }        
        
        
        timepoint_info = {};
        timepoint_info["div"] = {1,5}["0"];
        timepoint_info["pheno"] = {1,3}["0"];
        div_data = {1,5}["0"];
    }
    
        
    fprintf     (output_file, thisDate, "\t", Join("\t", timepoint_info["div"]), "\t", Join ("\t", timepoint_info["pheno"]), 
        "\t", Join ("\t", div_data) , "\n");
            
    ExecuteCommands ("DeleteObject (sampled`thisDate`_lf);");
}

//----------------------------------------------------------------------------------------

function fitBranchLength (s, ns) {
    //fprintf (stdout, tn, bn, s, ns, "\n");
    ExecuteCommands ("`tn`.`bn`.synRate = s");
    ExecuteCommands ("`tn`.`bn`.nonSynRate = ns");
    return Eval ("BranchLength (`tn`, \"`bn`\")");
}

function makePartitionBuiltTreeFitMG (dataID, sites, sequences, prefix, reroot) {
    DataSetFilter filteredData = CreateFilter(*dataID,1,sites,sequences);
     _normalizeSequenceNamesInAFilter ("filteredData");
    Export (fd, filteredData);
    DataSet reduced = ReadFromString (fd);
    
    UseModel (MGLocal);
    ExecuteCommands ("Tree `prefix`_tree = " + InferTreeTopology (1) + ";
    if (reroot != None) {
        tree_s = RerootTree (`prefix`_tree, reroot);
        Tree `prefix`_tree  = tree_s;
    }
    
    /*bls = getNucRevBranchLengthsAndParameters (\"filteredData\", \"`prefix`_tree\");
    
    bnames = BranchName (`prefix`_tree, -1);
    for (b = 0; b < Columns (bnames) - 1; b += 1) {
        tn = \"`prefix`_tree\";
        bn = bnames[b];
        v  = (bls[\"lengths\"])[bn];
        FindRoot (z, fitBranchLength (x, x*0.25) - v, x, 0, 10000);
    }
    
    vars = Rows (bls);
    for (k = 0; k < Columns (vars); k+=1) {
        if (vars[k] != \"lengths\") {
            ExecuteCommands (vars[k] + \"=\" + bls[vars[k]]);
        }
    }
    USE_LAST_RESULTS = 1;
    */
    
    DataSetFilter `prefix`_filter = CreateFilter(reduced,3,,,GeneticCodeExclusions);
    LikelihoodFunction `prefix`_lf = (`prefix`_filter,`prefix`_tree);
    unconstrainGlobalParameters(\"`prefix`_lf\");
    VERBOSITY_LEVEL = 0;
    Optimize (res, `prefix`_lf);    
    VERBOSITY_LEVEL = 0;
    USE_LAST_RESULTS = 0;
    ");  
    
    if (reroot != None) {
        if (filteredData.species > 2) {
            ExecuteCommands ("DataSet mrca_seqs = ReconstructAncestors (`prefix`_lf);");
            DataSetFilter mrca_filter = CreateFilter (mrca_seqs, 1);
            GetDataInfo (rootSeq, mrca_filter, 0);
        } else {
            GetDataInfo (rootSeq, filteredData, 0);            
        }    
        return rootSeq;
    }
    
    fixGlobalParameters ("`prefix`_lf");
    dSdN = _computeSNSSites ("`prefix`_filter", _Genetic_Code, vectorOfFrequencies, call_count);
    dS = dSdN["Sites"]/dSdN ["SSites"];
    dN = dSdN["Sites"]/dSdN ["NSSites"];
    return {"div" : treeDiv ("`prefix`_tree", dS, dN), "pheno": phenotypeAFilter ("`prefix`_filter"), "mrca": mrca};
    
}

//----------------------------------------------------------------------------------------

function treeDiv (treeiD,dS,dN) {
    leafCount = computeMultFactorsWeighted  (treeiD);
	BRANCH_LENGTH_STENCIL = _blStencils["Syn"];
	divInfo 		=	 computeTotalDivergence (treeiD);
	syn  			= 	2*divInfo[0]/leafCount/(leafCount-1);
	BRANCH_LENGTH_STENCIL = _blStencils["NonSyn"];
	divInfo 		=	 computeTotalDivergence (treeiD);
	ns  			= 	2*divInfo[0]/leafCount/(leafCount-1);
    BRANCH_LENGTH_STENCIL = 0;
    return {{syn__,ns__,syn__+ns__,syn__*dS__,ns__*dN}};    
}

//----------------------------------------------------------------------------------------

function divergenceForFilter (filterID, rootSeq) {
    total = {{0,0,0,0,0}};
    /*
    filter_root = *"filteredData.site_map";
    abridged_rootSeq = ""; abridged_rootSeq * 128;
    for (k = 0; k < Columns(filter_root); k+=1) {
        abridged_rootSeq * rootSeq[filter_root[k]];
    }
    abridged_rootSeq * 0;
    */
    
    copy_count = 0;
    
    for (k = 0; k < *"`filterID`.species"; k+=1) {
        GetDataInfo (seq, *filterID, k);
        GetString (seqName, *filterID, k);
        copyN = get_copy_number (seqName);
        pc = pairwiseCodon (rootSeq, seq);
        total = total + (pc)*(0+copyN);
        copy_count += copyN;
        //fprintf (stdout, ">",seqName, "\n", seq, "\n>root\n", rootSeq, "\n", pc[0], "\n");
        
        if (Abs (div_filter_name)) {
            div_filter_name_l = div_filter_name + "_" + seqName + ".lf";
            fprintf (div_filter_name_l, CLEAR_FILE, pairFunction);
        }
    }
    return total*(1/copy_count);
}

//----------------------------------------------------------------------------------------

function pairwiseCodon (seq1, seq2) {
    filter_string = "\n>1\n`seq1`\n>2\n`seq2`";
    DataSet pair = ReadFromString (filter_string);
    DataSetFilter pairFilter = CreateFilter (pair, 3, "", "", GeneticCodeExclusions);
    Tree pairTree = (1,2);
    LikelihoodFunction pairFunction = (pairFilter, pairTree);
    Optimize (res, pairFunction);
    dSdN = _computeSNSSites ("pairFilter", _Genetic_Code, vectorOfFrequencies, call_count);
    dS = dSdN["Sites"]/dSdN ["SSites"];
    dN = dSdN["Sites"]/dSdN ["NSSites"];
    return treeDiv ("pairTree",dS,dN);
}

//----------------------------------------------------------------------------------------

lfunction phenotypeASequence (seq) {
    //fprintf (stdout, "phenotypeASequence\n");

    l = Abs (seq^{{"[\\?X\\-]",""}});
    pngs = countPNGS (seq);
    iep = isoElectricPoint (seq);
    
    return {"Length": l, "PNGS": pngs, "Isoelectric Point": iep};
}


//----------------------------------------------------------------------------------------

lfunction phenotypeAFilter (filterName) {
    codon_2_aa = defineCodonToAA ();
    
    upper_limit = ^"`filterName`.species";
    all_phenotypes  = {upper_limit, 3};
    counts = {upper_limit,1};
    mean_phenotypes = {1, 3};
    
    
    for (seq_id = 0; seq_id < upper_limit; seq_id += 1) {
        GetDataInfo (this_sequence, ^filterName,seq_id);
        GetString   (this_seq_name, ^filterName,seq_id);
         
        aa_seq = translateCodonToAA(this_sequence,codon_2_aa,0);
        region_pheno             = phenotypeASequence(aa_seq);
        all_phenotypes [seq_id][0] = region_pheno["Length"];
        all_phenotypes [seq_id][1] = region_pheno["PNGS"];
        all_phenotypes [seq_id][2] = region_pheno["Isoelectric Point"];
        counts [seq_id] = get_copy_number (this_seq_name);
    
    }
    
    total = +counts;
    for (k = 0; k < 3; k+=1) {
        mean_phenotypes [k] = (+((all_phenotypes[-1][k])$counts))/total;
    }

    return mean_phenotypes;
}

