MESSAGE_LOGGING = 0;

LoadFunctionLibrary         ("chooseGeneticCode", {"0" : "Universal"});
LoadFunctionLibrary         ("ReadDelimitedFiles");
LoadFunctionLibrary         ("NJ");
LoadFunctionLibrary         ("BranchLengthFitters");
LoadFunctionLibrary         ("TreeTools");
LoadFunctionLibrary         ("LocalMGREV");
LoadFunctionLibrary         ("CF3x4");
LoadFunctionLibrary         ("TreeFunctions");
LoadFunctionLibrary         ("dSdNTreeTools");
LoadFunctionLibrary         ("_cwf.ibf");


fscanf                      (stdin, "String", filename);

DataSet ds                  = ReadDataFile (filename);

if (ds.species > 1) { 
    DataSetFilter filteredData          = CreateFilter (ds, 3, "", "", GeneticCodeExclusions);

    invmap = _normalizeSequenceNamesInAFilter ("filteredData");
    tree_string = InferTreeTopology (1);
    HarvestFrequencies (nuc3x4, filteredData, 3, 1, 1);
    nucCF						= CF3x4	(nuc3x4, GeneticCodeExclusions);
    ModelMatrixDimension        = CountSenseCodons (_Genetic_Code);
    codon3x4					= BuildCodonFrequencies (nucCF);
    PopulateModelMatrix			  ("MGMatrix",  nucCF);
    Model		MGlocal			= (MGMatrix, codon3x4, 0);
    Tree T = tree_string;

    LikelihoodFunction lf = (filteredData, T);
    Optimize (res, lf);

    theStencil = ComputeScalingStencils (0);

    result = {"Diversity": {}, "BranchLength": {}, "Tree": {}};

    BRANCH_LENGTH_STENCIL = theStencil["Syn"];
    sR = getCI ("T", 0);
    (result["Diversity"])["S"] = sR[0];
    (result["BranchLength"])["S"] = sR[1];
    (result["Tree"])["S"] = Format (T,1,1);

    BRANCH_LENGTH_STENCIL = theStencil["NonSyn"];
    nsR = getCI ("T", 0);

    (result["Diversity"])["NS"]    = nsR[0];
    (result["BranchLength"])["NS"] = nsR[1];
    (result["Tree"])["NS"] = Format (T,1,1);


    BRANCH_LENGTH_STENCIL = 0;
    (result ["Tree"])["Combined"] = Format (T, 1, 1);

    global R = 1;
    ReplicateConstraint ("this1.?.nonSynRate:=R*this2.?.synRate", T,T);
    USE_LAST_RESULTS = 1;
    Optimize (res, lf);

    result ["omega"] = R;
} else {
   result = {"Diversity": {"S": 0, "NS" : 0}, "BranchLength": {"S": 0, "NS" : 0}, "Tree": "", "omega": 0};
}

fprintf (stdout, result, "\n");



/*----------------------------------------------------------------*/

function getCI (treeName, doCI) {
	
	leafCount = computeMultFactorsWeighted  (treeName);


	divInfo 		=	 computeTotalDivergence (treeName);
	pInfo 			= 	2*divInfo[0]/leafCount/(leafCount-1);
	currentDepth	= 	divInfo[1]/(Abs(treeAVL2)-2);
		
	if (!doCI) {
	    return {{pInfo, currentDepth}};
	}

	lf_Count						    = Rows ("LikelihoodFunction");
    for (lf_ID = 0; lf_ID < lf_Count; lf_ID = lf_ID + 1)
    {
        GetString (treeID, LikelihoodFunction,lf_ID);
        GetString(lfInfo,^treeID,-1);
        lfTrees = lfInfo["Trees"];
        for (k = 0; k<Columns(lfTrees); k=k+1)
        {
            if (lfTrees[k] == treeName)
            {
                break;
            }
        }
        if (k < Columns(lfTrees))
        {
            //fprintf (stdout, "\nTree ",Columns(lfTrees)," is a part of likelihood function ", treeID, "\n");
            global			TreeScalerParameter = 1;
            COVARIANCE_PARAMETER = "TreeScalerParameter";
            ExecuteCommands	("ClearConstraints("+treeName+");ReplicateConstraint (\"this1.?.?:=TreeScalerParameter*this2.?.?__\","+treeName+","+treeName+");\n");
            
            COVARIANCE_PRECISION = 0.95;
            ExecuteCommands ("CovarianceMatrix(cmx,"+treeID+");");
            fprintf (stdout, "\nMultiplicative range ", cmx[0], "-", cmx[2], "\n");
            fprintf (stdout, "Mean pairwise divergence: ", pInfo, " (", pInfo*cmx[0],",",pInfo*cmx[2],")\n");
            fprintf (stdout, "Mean branch length divergence: ", currentDepth, " (", currentDepth*cmx[0],",",currentDepth*cmx[2],")\n");
            
            return {"MEAN PAIRWISE DIVERGENCE": {{pInfo, pInfo*cmx[0], pInfo*cmx[2]}},
                    "MEAN BRANCH LENGTH": {{currentDepth, currentDepth*cmx[0],currentDepth*cmx[2]}}}; 
        }
    }
    return 0;
}