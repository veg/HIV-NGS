#!/usr/bin/env julia

#------------------------------------------------------------------------------#
#                                 mcmc.jl                                      #
#------------------------------------------------------------------------------#

# First, import some dependencies.

using ArgParse
using FastaIO
using JSON

# Assign the AMBIGS dictionary.

const AMBIGS = Dict{String, Char}(
    [
        "A", "C", "G", "T",
        "AC", "AG", "AT", "CG", "CT", "GT",
        "ACG", "ACT", "AGT", "CGT",
        "ACGT", ""],
    [
        'A', 'C', 'G', 'T',
        'M', 'R', 'W', 'S', 'Y', 'K',
        'V', 'H', 'D', 'B',
        'N', '-'])




#------------------------------------------------------------------------------#
#                          functional definitions                              #
#------------------------------------------------------------------------------#




#----------------------------------index---------------------------------------#

# Define a function that returns all indices at which a string matches a 
# certain character.

function index(s::String, c::Char)

    # For each character in the string, check if it matches the input character
    # and if it does, return the index.

    for i in 1:length(s)
        @inbounds s_i = s[i]
        if c == s_i
            return i
        end
    end
    
    # Return 0 if there are no matches.
    
    0
end




#---------------------------------countmsa-------------------------------------#

# Define a function to count, for each character of an alphabet, the number
# of times that character appears in each sequence of an MSA file.

function countmsa(msa::String, alphabet::String)

    # Initialize the counters, counts and ncols.
    
    counts = nothing
    ncols = 0
    tic()
    
    # For each sequence in the input .msa file,
    
    for (i, (name, seq)) in enumerate(FastaReader(msa))
    
        # If the counters are empty, assign to ncols the length of the current
        # sequence and initialize counts as a vector of zeros with length ncols.
    
        if counts == nothing
            ncols = length(seq)
            counts = zeros(Int64, (ncols, 4))
            
        # If the counters are not empty, check that ncols == length(seq), and
        # say something if they are not equal.
            
        elseif length(seq) != ncols
            error("provided file is not an MSA! length($name) = $(length(seq)), not $ncols!")
        end
        
        # For each index of the current sequence,
        
        for j in 1:length(seq)
        
           # Check which character is at the current index and count it.
        
            @inbounds L = uppercase(seq[j])
            for k in 1:length(alphabet)
                @inbounds c = alphabet[k]
                if L == c
                    counts[j, k] += 1
                end
            end
        end
        progress("loading msa: read $i sequences .. ")
    end
    done()
    
    # Return the array of counts.
    
    counts
end




#---------------------------------rategrid-------------------------------------#

# Define a function to generate a rate grid of a given size and error threshold.

function rategrid(n::Integer, error_threshold::Float64)
    update("generating rate grid .. ")
    
    # Generate a vector of rates increasing on a log scale.
    
    if n == 10
        rates = [
            0.0001, 0.0005, 0.001, 0.002, 0.004, 0.008 , 0.01, 0.02, 0.05, 0.25
        ]
    else   
        rates = logspace (-4,log10(0.1),n)
    end 
    #print (rates)
    
    # Initialize uniques, an empty two-dimensional array. Uniques is initialized
    # as a set to avoid duplicate rows.
    
    uniques = Set{Vector{Float64}}()
    
    # Build uniques using the error threshold.
    
    M = 1. - 3. * error_threshold
	union! (uniques, {[M,error_threshold,error_threshold,error_threshold]})
	union! (uniques, {[error_threshold,M,error_threshold,error_threshold]})
	union! (uniques, {[error_threshold,error_threshold,M,error_threshold]})
	union! (uniques, {[error_threshold,error_threshold,error_threshold,M]})
        
    # For each value in the rates vector,
        
    for a in rates
    
        # Add a row to the uniques array.
    
        M = 1. - a - 2. * error_threshold
        
	    for i1 in 1:4
	        for i2 in 1:4
	            if i1 != i2
                    rate_vec = fill (error_threshold, (4))
	                rate_vec[i1] = M
	                rate_vec[i2] = a
	                union! (uniques, {rate_vec})
	            end
	        end
	    end
	    
	    # Add four more rows to the uniques array.
                
        for b in rates
            for c in rates
                M = 1. - (a + b + c)
		        union! (uniques, {[M,a,b,c]})
		        union! (uniques, {[a,M,b,c]})
                union! (uniques, {[a,b,M,c]})
                union! (uniques, {[a,b,c,M]})

            	
            end
        end
    end
    done()
    
    # Initialize a counter i.
    
    i = 1
    
    # Initialize an empty rate grid.
    
    rg = Array (Float64, (length (uniques),4))
    
    # For each row in uniques,
    
    for r in uniques
	    #println (r)
	    
	    # Set the present row of the rate grid equal to the present row of
	    # uniques.
	    
	    rg[i, :] = r 
	    
	    # Advance the counter i.
	    
	    i += 1
    end
    
    # Print the length of uniques and return the rate grid.
    
    println (length (uniques))
    #print (uniques)
    return rg
end




#----------------------------------lfact---------------------------------------#

# Initialize the lfact_cache

lfact_cache = Float64[0.0]

# Define a function to return the lfact (logarithm of factorial) for a given
# integer.

function lfact(N::Integer)

    # If the input integer is 0, return 0.

    if N == 0
        return 0.0
    end
    
    # Otherwise, note the current length of the lfact_cache and check if it is
    # less than the input integer.
    
    n = length(lfact_cache)
    if n < N
    
        # If so, extend the lfact_cache to be of length equal to the input 
        # integer and fill the new space with logarithmically-growing values.
    
        resize!(lfact_cache, N)
        for i in (n + 1):N
            @inbounds lfact_cache[i] = lfact_cache[i - 1] + log(i)
        end
    end
    
    # Return the largest value of the lfact_cache.
    
    @inbounds return lfact_cache[N]
end


#-----------------------------------lmc----------------------------------------#


# Define a function to return the lmc value for a given site.

function lmc(counts::Array{Int64,2}, site::Int64, nchars::Int64)

    # Initialize placeholders.

    c = 0.0
    s = 0

    # For each character,
    
    for char in 1:nchars
        @inbounds c -= lfact(counts[site, char])
        @inbounds s += counts[site, char]
    end
    return c + lfact(s)
end




#--------------------------------gridscores------------------------------------#

# Approximate negative infinity.

const smin = log(realmin(Float64))

# Define a function to compute grid scores.

function gridscores(counts::Array{Int64,2}, rates::Array{Float64,2}, alphabet::String)
    update("computing grid scores .. ")
    nchars = length(alphabet)
    npoints, _ = size(rates)
    nsites, _ = size(counts)
    
    # Initialize arrays for conditionals and scalers.
    
    conditionals = Array(Float64, (npoints, nsites))
    scalers = Array(Float64, (nsites,))
    
    # Compute log (rates)
    
    log_rates = log (rates)
     
    # For each site of the input .msa, 
    
    for i in 1:nsites
    
        # Set m to approximately negative infinity
    
        m = -realmax(Float64)
        
        # Compute lmc for the current site
        
        @inbounds c = lmc(counts, i, nchars)
        
        # For each rate,
        
        for j in 1:npoints
           s = c
           
           # for each character,
           
           for k in 1:nchars
                @inbounds s += counts[i, k] * log_rates[j, k]
            end
            if s > m
                m = s
            end
            
            # Update the conditionals.
            
            @inbounds conditionals[j, i] = s
        end
        
        # Update the scalers.
        
        @inbounds scalers[i] = m
        
        # For each rate,
        
        for j in 1:npoints
            @inbounds s = conditionals[j, i] - m
            if s < smin
                @inbounds conditionals[j, i] = 0
            else
                @inbounds conditionals[j, i] = exp(s)
            end
            #println (conditionals [j, i])
            #if (i == 771) && conditionals[j, i] > 0
            #    println (conditionals[j, i], "\t", join (rates[j , :], "\t"))
            #end
        end
        
            
    end
    done()
    
    # Return the conditionals and scalers.
    
    return (conditionals, scalers)
end




#---------------------------------ldensity-------------------------------------#

# Define a function to compute ldensity.

function ldensity(weights::Array{Float64,2}, alpha::Float64)

    # Set bounds on the weights array and on alpha to protect against extreme
    # cases.

    if minimum(weights) <= 1e-10
        return 1e-10
    end
    if alpha == 1
        return 0.0
    end
    
    # Compute the density.
    
    n = length(weights)
    r = lgamma(n * alpha) - n * lgamma(alpha)
    for i in 1:n
        @inbounds r += log(weights[i]) * (alpha - 1)
    end
    
    # Return the value of the density given weights and alpha.
    
    return r
end




#------------------------------------jll---------------------------------------#

# Define a function to compute ll and dir.

function jll(conditionals::Array{Float64,2}, scalers::Array{Float64,1}, weights::Array{Float64,2}, alpha::Float64)
    ll = sum(log(weights * conditionals) + scalers')
    dir = ldensity(weights, alpha)
    return (ll, dir)
end




#-----------------------------------mcmc---------------------------------------#

# Define a function to carry out the mcmc process. Alpha was once called
# concentration_parameter.

function mcmc(
        conditionals::Array{Float64,2},
        scalers::Array{Float64,1};
        chain_length::Int64=10_000_000,
        burnin_fraction::Float64=0.5,
        expected_nsamples::Int64=100,
        alpha::Float64=0.5)
        
    # Compute nburnin from the inputs.

    nburnin = iround(chain_length * burnin_fraction)
    
    # Note the shape of conditionals.
    
    npoints, nsites = size(conditionals)
    
    # Normalize the conditionals and weights by site.

    normalized_by_site = ones((1, npoints)) * conditionals
    
    normalized_weights = conditionals * ((1 ./ normalized_by_site) .* eye(nsites))
    
    # Sum the weights by site and normalize.
    
    sum_by_site = normalized_weights * ones((nsites, 1))
    weights = transpose((rand((npoints, 1)) .* (1.2 * sum_by_site - 0.8 * sum_by_site)) + 0.8 * sum_by_site)
    weights /= sum(weights)
    
    # Determine step size.

    stepsize = max (2. / max(npoints, nsites), median(weights))
    
    # Determine the number of samples and initialize arrays for their weights
    # and likelihoods.

    nsample = div(chain_length - nburnin, expected_nsamples)
    sampled_weights = zeros(Float64, (expected_nsamples, npoints))
    sampled_likelihoods = zeros(Float64, (expected_nsamples,))
    
    # Compute likelihoods, logsum, and jll for the current site.

    current_site_likelihoods = weights * conditionals
    current_site_logsum = sum(log(current_site_likelihoods))
    current_ll, current_dir = jll(conditionals, scalers, weights, alpha)
    
    # Initialize step trackers.

    step = 1
    accepted_steps = 0
    mean_sampled_ll = 0
    sample_idx = 0
    
    # Initialize some other objects.

    diffvec = Array(Float64, (1, nsites))
    idxs = 1:npoints
    llstr = "(burning in)"

    tic()
    
    # For each step of the mcmc simulation,
    
    for step in 1:chain_length

        # Choose two random but separate indices.

        i = rand(idxs)
        j = rand(idxs)
        while i == j
            j = rand(idxs)
        end
        
        # Choose a random number proportional to stepsize and choose two
        # of the weights.

        change = rand() * stepsize
        @inbounds weights_i = weights[i]
        @inbounds weights_j = weights[j]
        
        # If weights_i exceeds change,

        if weights_i > change
        
            # Initialize lldiff
            
            lldiff = 0
            
            # For each site, 
            
            for k in 1:nsites
            
                # Add a corresponding entry to diffvec and update lldiff.
            
                @inbounds diffvec[k] = (conditionals[j, k] - conditionals[i, k]) * change
                @inbounds lldiff += log(current_site_likelihoods[k] + diffvec[k])
            end
            
            # Normalize lldiff, generate dirdiff, and infer movecost.
            
            lldiff -= current_site_logsum
            @inbounds dirdiff = (alpha - 1) * (log((weights_i - change) / weights_i) + log((weights_j + change) / weights_j))
            movecost = exp(lldiff + dirdiff)
            
            # Use movecost to determine the probability of changing state.

            if rand() <= movecost
            
                # Update values and weights for the current site.
            
                current_ll += lldiff
                current_dir += dirdiff

                current_site_likelihoods += diffvec
                current_site_logsum += lldiff

                @inbounds weights[i] -= change
                @inbounds weights[j] += change
                
                # Advance the accepted_steps counter.

                accepted_steps += 1
            end
        end
        
        # If this is not a burn-in step, 

        if step > nburnin
        
            # If a sampling interval has elapsed,
            
            if (step - nburnin) % nsample == 0
            
                # Advance the sample index counter, and sample the weights and
                # likelihoods from this simulation step.
            
                sample_idx += 1
                sampled_weights[sample_idx, :] = weights
                sampled_likelihoods[sample_idx] = current_ll
                mean_sampled_ll += current_ll
                llstr = "$(trunc(mean_sampled_ll / sample_idx, 3))"
            end
        end
        
        # If a sampling interval has elapsed,

        if step % nsample == 0
        
            # Report progress.
        
            progress("running MCMC chain: step = $step/$chain_length, mean logL = $llstr, acceptance rate = $(trunc(accepted_steps / step, 3)) .. ")
        end
    end

    done()
    
    # Return all of the sampled likelihoods and weights.

    return (sampled_likelihoods, sampled_weights)
end




#---------------------------------loadrates------------------------------------#

# Define a function to load the rates file.

function loadrates()
    rates3 = readdlm("data/rates.txt", ' ', Float64)
    nrates, _ = size(rates3)
    rates = Array(Float64, (nrates, 4))
    rates[:, 1:3] = rates3
    for i in 1:nrates
        @inbounds rates[i, 4] = 1.0 - sum(rates[i, 1:3])
    end
    return rates
end




#---------------------------------postproc-------------------------------------#

# Define a function to compute prior and posterior probabilities for the 
# conditionals and ws arrays.

function postproc(conditionals::Array{Float64,2}, ws::Array{Float64,2})
    update("computing posterior probabilities .. ")
    
    # Record the shapes of conditionals and ws.
    
    nsites, npoints = size(conditionals)
    nsamples, _ = size(ws)
    
    # Initialize an array to store the priors.
    
    priors = zeros(Float64, (1, npoints))
    
    # Deduce the priors from ws and normalize.
    
    for i in 1:nsamples
        @inbounds priors += ws[i, :]
    end
    priors = priors' / nsamples
    priors /= sum(priors)
    
    # Compute the normalization constant for the posteriors.
    
    normalization = conditionals * priors
    
    # Compute the posteriors.
    
    posteriors = ((1 ./ normalization) .* eye(nsites)) * (conditionals * (priors .* eye(npoints)))
    done()
    
    # Return the priors and posteriors.
    
    return (priors, posteriors)
end




#-------------------------------callvariants-----------------------------------#

# Define a function to make variant calls.

function callvariants(
        counts::Array{Int64,2},
        rates::Array{Float64,2},
        priors::Array{Float64,2},
        posteriors::Array{Float64,2},
        alphabet::String;
        target_rate::Float64=0.1,
        posterior_threshold::Float64=0.999,
        json_file=nothing,
        min_coverage::Int64=100)

    update("calling variants .. ")
    
    # Record the number of posteriors and filter them according to the target
    # rate.
    
    nsites, npoints = size(posteriors)
    posterior_prob = posteriors * map(x -> x > target_rate, rates)
    majority_prob  = posteriors * map(x -> x >= 0.5, rates)
    
    # If json_file already exists, 
    
    if json_file != nothing
    
        # Initialize the overall report and rate dictionaries.
        
        overall_report = Dict()
        rate_dict = Dict {Int64, Dict{String,Float64}}()
        
        # For each posterior,
        
        for i in 1:npoints
        
            # If the corresponding prior is great enough,
        
            if priors[i] >= 0.5/nsites
            
                # Add the rates and prior at this index to the rate dictionary.
            
                rate_dict [i] = (String=>Float64)[
                        string(alphabet[N]) => rates[i,N] for N = 1:length(alphabet) ]
                rate_dict [i]["weight"] = priors[i]
            end 
        end
        
        # Record some more information.
        
        posterior_means = posteriors * rates
        posterior_dict = Dict {Int64, Dict{String,Any}}()
        overall_report ["priors"] = rate_dict
        consensus = ""
    end
    
    # Initialize the variants dictionary.
    
    variants = Dict{Int64,String}()
    
    # For each site,
    
    for i in 1:nsites
        s = Set{Char}()
        
        # Record the majority pair.
        
        maj_pair = findmax (majority_prob[i, :])
        
        # If the majority posterior is large enough,
        
        if maj_pair[1] >= posterior_threshold
        
            # Add the majority prior to s.
            
            push!(s, alphabet[maj_pair[2]])
            
            # For each character in the alphabet,
            
            for j in 1:length(alphabet)
            
                # If the posterior at the current character is large enough,
            
                if posterior_prob[i, j] >= posterior_threshold
                
                    # Add the current character to s.
                
                    push!(s, alphabet[j])
                end
            end
        end
        
        # Organize s, and assign it to the variants at this index.
        
        s_ = join(sort!(collect(s)))
        #if length(s) > 1
        println("$i: $s_,\t[$(join(counts[i, :], '\t'))]")
        #end
        variants[i] = s_
        
        # If json_file already exists,
        
        if json_file != nothing
        
            # If coverage at this index is above the minimum threshold,
        
            if sum (counts[i,:]) >= min_coverage
            
                # Record the index of the maximum posterior probability and 
                # update consensus.
            
                index = indmax (posterior_prob[i,:])
                consensus = *(consensus,string(alphabet[index]))
            else
            
                # Apply a null update to consensus.
            
                consensus = *(consensus,"-")
            end
            
            # Update the posterior dictionary.
            
            posterior_dict[i] = (String=>Any)[
                string(alphabet[N]) => [posterior_prob[i,N], posterior_means[i,N], counts[i,N]] for N = 1:length(alphabet) ]
            posterior_dict[i]["variants"] = s_
        end 
    end
    done()
    
    # If json_file already exists,    

    if json_file != nothing
    
        # update the overall report and write it to json_file.
    
        overall_report["consensus"] = consensus
        overall_report["posteriors"] = posterior_dict
        json_out = open (json_file, "w")
        JSON.print (json_out, overall_report)
        close (json_out)
    end
    
    # Return the variants.
    
    return variants
end




#--------------------------------filterseqs------------------------------------#

# Define a function to filter the sequences.

function filterseqs(msa::String, variants::Dict{Int64,String}, dest::String, alphabet::String)
    tic()
    
    # Extract a sites list from the variants.
    
    sites = sort!(collect(keys(variants)))
    
    # Initialize varmask.
    
    varmask = zeros(Bool, ((length(variants), 4)))
    
    # For each site,
    
    for i in 1:length(sites)
        site = sites[i]
        
        # For each character in the alphabet,
        
        for j in 1:length(alphabet)
        
            # If the current character is a variant at the current site,
        
            if alphabet[j] in variants[site]
            
               # Record it in varmask.
            
                varmask[i, j] = true
            end
        end
    end
    
    # Open a fasta file.
    
    FastaWriter(dest == "-" ? STDOUT : dest) do fw
    
        # For each sequence in the input .msa file,
    
        for (i, (name, seq)) in enumerate(FastaReader(msa))
        
            # Initialize the counter v,
        
            v = 1
            seq_ = Array(Char, (length(seq),))
            
            # For each index in the current sequence,
            
            for j in 1:length(seq)
                seq_[j] = seq[j]
                
                # If the sequence has not ended, and the counter v is equal
                # to the current indes,
                
                if v < length(sites) && j == sites[v]
                
                    # Extract the alphabetical index of the character at the
                    # current index.
                
                    idx = index(alphabet, seq[j])
                    
                    # If the alphabetical index is positive and this character
                    # is not a variant at this site, but the variant at this 
                    # site is ambiguous,
                    
                    if idx > 0 && !varmask[v, idx] && variants[sites[v]] in keys(AMBIGS)
                    
                        # Assign the appropriate AMBIGS element to the current
                        # index.
                    
                        seq_[j] = AMBIGS[variants[sites[v]]]
                    end
                    
                    # Advance the counter v.
                    
                    v += 1
                end
            end
            
            # Report progress.
            
            writeentry(fw, name, join(seq_, ""))
            progress("filtering unsupported variants: processed $i sequences .. ")
        end
    end
    done()
end




#------------------------------------------------------------------------------#
#                         completion announcement                              #
#------------------------------------------------------------------------------#

# Report

function done()
    print_with_color(:blue, STDERR, "done, took $(toq())s")
    println(STDERR)
end
function progress(args...)
    print_with_color(:blue, STDERR, "\rINFO: ", args...)
end
function update(args...)
    print_with_color(:blue, STDERR, "INFO: ", args...)
    tic()
end




#------------------------------------------------------------------------------#
#                                   main                                       #
#------------------------------------------------------------------------------#

# Define the main function.

function main(args)

    # Parse arguments.

    s = ArgParseSettings()
    s.description = "call variants using a multinomial model sampled by MCMC"

    @add_arg_table s begin
        "--grid-density", "-g"
            arg_type = Int64
            default = 10
        "--chain-length", "-c"
            arg_type = Int64
            default = 2_000_000
        "--burnin-fraction", "-b"
            arg_type = Float64
            default = 0.5
        "--target-rate", "-t"
            arg_type = Float64
            default = 0.01
        "--posterior-threshold", "-p"
            arg_type = Float64
            default = 0.95
        "--filter", "-f"
            arg_type = String
        "--json-report", "-j"
            arg_type = String
        "msa"
            required = true
    end

    parsed_args = parse_args(args, s)
    
    # Construct the rate grid.
    
    rates = rategrid(parsed_args["grid-density"], parsed_args["target-rate"])
    
    # Construct the alphabet.

    const alphabet = "ACGT"
    
    # Execute countmsa() on the input .msa file and the alphabet.
        
    counts = countmsa(parsed_args["msa"], alphabet)
    
    # Construct conditionals and scalers using gridscores().

    conditionals, scalers = gridscores(counts, rates, alphabet)
    
    # Obtain lls and ws via MCMC simulation.

    lls, ws = mcmc(
        conditionals, scalers,
        chain_length=parsed_args["chain-length"], burnin_fraction=parsed_args["burnin-fraction"])
        
    # Extract priors and posteriors from conditionals and ws using postproc.

    priors, posteriors = postproc(conditionals', ws)
    
    # Call variants.

    variants = callvariants(
        counts, rates, priors, posteriors, alphabet,
        target_rate=parsed_args["target-rate"], posterior_threshold=parsed_args["posterior-threshold"], json_file = parsed_args["json-report"])
    
    # If instructed to, filter the sequences.
    
    if "filter" in keys(parsed_args)
        filterseqs(parsed_args["msa"], variants, parsed_args["filter"], alphabet)
    end
end

# Execute.

main(ARGS)
