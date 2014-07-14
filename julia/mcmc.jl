#!/usr/bin/env julia

using ArgParse
using FastaIO
using JSON

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

function index(s::String, c::Char)
    for i in 1:length(s)
        @inbounds s_i = s[i]
        if c == s_i
            return i
        end
    end
    0
end

function countmsa(msa::String, alphabet::String)
    counts = nothing
    ncols = 0
    tic()
    for (i, (name, seq)) in enumerate(FastaReader(msa))
        if counts == nothing
            ncols = length(seq)
            counts = zeros(Int64, (ncols, 4))
        elseif length(seq) != ncols
            error("provided file is not an MSA! length($name) = $(length(seq)), not $ncols!")
        end
        for j in 1:length(seq)
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
    counts
end

function rategrid(n::Integer, error_threshold::Float64)
    update("generating rate grid .. ")
    if n == 10
        rates = [0.0001, 0.0005, 0.001, 0.002, 0.004, 0.008 , 0.01, 0.02, 0.05, 0.25]
    else   
        rates = logspace (-4,log10(0.1),n)
    
    end 
    print (rates)
    uniques = Set{Array{Float64, 2}}()
    
    M = 1. - 3. * error_threshold
	add! (uniques, [M error_threshold error_threshold error_threshold])
	add! (uniques, [error_threshold M error_threshold error_threshold])
	add! (uniques, [error_threshold error_threshold M error_threshold])
	add! (uniques, [error_threshold error_threshold error_threshold M])
        
        
    for a in rates
        M = 1. - a - 2. * error_threshold
        
	    for i1 in 1:4
	        for i2 in 1:4
	            if i1 != i2
                    rate_vec = fill (error_threshold, (1,4))
	                rate_vec[i1] = M
	                rate_vec[i2] = a
	                add! (uniques, rate_vec)
	            end
	        end
	    end
                
        for b in rates
            for c in rates
                M = 1. - (a + b + c)
		        add! (uniques, [M a b c])
		        add! (uniques, [a M b c])
                add! (uniques, [a b M c])
                add! (uniques, [a b c M])
            	
            end
        end
    end
    done()
    i = 1
    
    
    rg = Array (Float64, (length (uniques),4))
    for r in uniques
	    #print (r)
	    rg[i, :] = r 
	    i += 1
    end
    println (length (uniques))
    return rg
end

lfact_cache = Float64[0.0]

function lfact(N::Integer)
    if N == 0
        return 0.0
    end
    n = length(lfact_cache)
    if n < N
        resize!(lfact_cache, N)
        for i in (n + 1):N
            @inbounds lfact_cache[i] = lfact_cache[i - 1] + log(i)
        end
    end
    @inbounds return lfact_cache[N]
end

function lmc(counts::Array{Int64,2}, site::Int64, nchars::Int64)
    c = 0.0
    s = 0
    for char in 1:nchars
        @inbounds c -= lfact(counts[site, char])
        @inbounds s += counts[site, char]
    end
    return c + lfact(s)
end

const smin = log(realmin(Float64))

function gridscores(counts::Array{Int64,2}, rates::Array{Float64,2}, alphabet::String)
    update("computing grid scores .. ")
    nchars = length(alphabet)
    npoints, _ = size(rates)
    nsites, _ = size(counts)
    conditionals = Array(Float64, (npoints, nsites))
    scalers = Array(Float64, (nsites,))
    log_rates = log (rates)
    for i in 1:nsites
        m = -realmax(Float64)
        @inbounds c = lmc(counts, i, nchars)
        for j in 1:npoints
           s = c
           for k in 1:nchars
                @inbounds s += counts[i, k] * log_rates[j, k]
            end
            if s > m
                m = s
            end
            @inbounds conditionals[j, i] = s
        end
        @inbounds scalers[i] = m
        for j in 1:npoints
            @inbounds s = conditionals[j, i] - m
            if s < smin
                @inbounds conditionals[j, i] = 0
            else
                @inbounds conditionals[j, i] = exp(s)
            end
            #if (i == 771) && conditionals[j, i] > 0
            #    println (conditionals[j, i], "\t", join (rates[j , :], "\t"))
            #end
        end
        
            
    end
    done()
    return (conditionals, scalers)
end

function ldensity(weights::Array{Float64,2}, alpha::Float64)
    if min(weights) <= 1e-10
        return 1e-10
    end
    if alpha == 1
        return 0.0
    end
    n = length(weights)
    r = lgamma(n * alpha) - n * lgamma(alpha)
    for i in 1:n
        @inbounds r += log(weights[i]) * (alpha - 1)
    end
    return r
end

function jll(conditionals::Array{Float64,2}, scalers::Array{Float64,1}, weights::Array{Float64,2}, alpha::Float64)
    ll = sum(log(weights * conditionals) + scalers')
    dir = ldensity(weights, alpha)
    return (ll, dir)
end

# alpha was once called concentration_parameter
function mcmc(
        conditionals::Array{Float64,2},
        scalers::Array{Float64,1};
        chain_length::Int64=10_000_000,
        burnin_fraction::Float64=0.5,
        expected_nsamples::Int64=100,
        alpha::Float64=0.5)

    nburnin = iround(chain_length * burnin_fraction)

    npoints, nsites = size(conditionals)

    normalized_by_site = ones((1, npoints)) * conditionals
    normalized_weights = conditionals * ((1 / normalized_by_site) .* eye(nsites))
    sum_by_site = normalized_weights * ones((nsites, 1))

    weights = transpose((rand((npoints, 1)) .* (1.2 * sum_by_site - 0.8 * sum_by_site)) + 0.8 * sum_by_site)
    weights /= sum(weights)

    stepsize = max (2. / max(npoints, nsites), median(weights))

    nsample = div(chain_length - nburnin, expected_nsamples)
    sampled_weights = zeros(Float64, (expected_nsamples, npoints))
    sampled_likelihoods = zeros(Float64, (expected_nsamples,))

    current_site_likelihoods = weights * conditionals
    current_site_logsum = sum(log(current_site_likelihoods))
    current_ll, current_dir = jll(conditionals, scalers, weights, alpha)

    step = 1
    accepted_steps = 0
    mean_sampled_ll = 0
    sample_idx = 0

    diffvec = Array(Float64, (1, nsites))
    idxs = 1:npoints
    llstr = "(burning in)"

    tic()

    for step in 1:chain_length

        i = rand(idxs)
        j = rand(idxs)
        while i == j
            j = rand(idxs)
        end

        change = rand() * stepsize
        @inbounds weights_i = weights[i]
        @inbounds weights_j = weights[j]

        if weights_i > change
            lldiff = 0
            for k in 1:nsites
                @inbounds diffvec[k] = (conditionals[j, k] - conditionals[i, k]) * change
                @inbounds lldiff += log(current_site_likelihoods[k] + diffvec[k])
            end
            lldiff -= current_site_logsum
            @inbounds dirdiff = (alpha - 1) * (log((weights_i - change) / weights_i) + log((weights_j + change) / weights_j))
            movecost = exp(lldiff + dirdiff)

            if rand() <= movecost
                current_ll += lldiff
                current_dir += dirdiff

                current_site_likelihoods += diffvec
                current_site_logsum += lldiff

                @inbounds weights[i] -= change
                @inbounds weights[j] += change

                accepted_steps += 1
            end
        end

        if step > nburnin
            if (step - nburnin) % nsample == 0
                sample_idx += 1
                sampled_weights[sample_idx, :] = weights
                sampled_likelihoods[sample_idx] = current_ll
                mean_sampled_ll += current_ll
                llstr = "$(trunc(mean_sampled_ll / sample_idx, 3))"
            end
        end

        if step % nsample == 0
            progress("running MCMC chain: step = $step/$chain_length, mean logL = $llstr, acceptance rate = $(trunc(accepted_steps / step, 3)) .. ")
        end
    end

    done()

    return (sampled_likelihoods, sampled_weights)
end

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

function postproc(conditionals::Array{Float64,2}, ws::Array{Float64,2})
    update("computing posterior probabilities .. ")
    nsites, npoints = size(conditionals)
    nsamples, _ = size(ws)
    priors = zeros(Float64, (1, npoints))
    for i in 1:nsamples
        @inbounds priors += ws[i, :]
    end
    priors = priors' / nsamples
    priors /= sum(priors)
    normalization = conditionals * priors
    posteriors = ((1 / normalization) .* eye(nsites)) * (conditionals * (priors .* eye(npoints)))
    done()
    return (priors, posteriors)
end

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
    
    
    nsites, npoints = size(posteriors)
    posterior_prob = posteriors * map(x -> x > target_rate, rates)
    majority_prob  = posteriors * map(x -> x >= 0.5, rates)
    
    if json_file != nothing
        overall_report = Dict()
        rate_dict = Dict {Int64, Dict{String,Float64}}()
        for i in 1:npoints
            if priors[i] >= 0.5/nsites
                rate_dict [i] = (String=>Float64)[
                        string(alphabet[N]) => rates[i,N] for N = 1:length(alphabet) ]
                rate_dict [i]["weight"] = priors[i]
            end 
        end
        posterior_means = posteriors * rates
        posterior_dict = Dict {Int64, Dict{String,Any}}()
        overall_report ["priors"] = rate_dict
        consensus = ""
    end
    
    variants = Dict{Int64,String}()
    for i in 1:nsites
        s = Set{Char}()
        
        maj_pair = findmax (majority_prob[i, :])
        if maj_pair[1] >= posterior_threshold
            push!(s, alphabet[maj_pair[2]])
            for j in 1:length(alphabet)
                if posterior_prob[i, j] >= posterior_threshold
                    push!(s, alphabet[j])
                end
            end
        end
        
        s_ = join(sort!(collect(s)))
        #if length(s) > 1
        println("$i: $s_,\t[$(join(counts[i, :], '\t'))]")
        #end
        variants[i] = s_
        if json_file != nothing
            if sum (counts[i,:]) >= min_coverage
                index = indmax (posterior_prob[i,:])
                consensus = *(consensus,string(alphabet[index]))
            else
                consensus = *(consensus,"-")
            end
            posterior_dict[i] = (String=>Any)[
                string(alphabet[N]) => [posterior_prob[i,N], posterior_means[i,N], counts[i,N]] for N = 1:length(alphabet) ]
            posterior_dict[i]["variants"] = s_
        end
        
    end
    done()
    
    if json_file != nothing
        overall_report["consensus"] = consensus
        overall_report["posteriors"] = posterior_dict
        json_out = open (json_file, "w")
        JSON.print (json_out, overall_report)
        close (json_out)
    end
    
    return variants
end

function filterseqs(msa::String, variants::Dict{Int64,String}, dest::String, alphabet::String)
    tic()
    sites = sort!(collect(keys(variants)))
    varmask = zeros(Bool, ((length(variants), 4)))
    for i in 1:length(sites)
        site = sites[i]
        for j in 1:length(alphabet)
            if alphabet[j] in variants[site]
                varmask[i, j] = true
            end
        end
    end
    FastaWriter(dest == "-" ? STDOUT : dest) do fw
        for (i, (name, seq)) in enumerate(FastaReader(msa))
            v = 1
            seq_ = Array(Char, (length(seq),))
            for j in 1:length(seq)
                seq_[j] = seq[j]
                if v < length(sites) && j == sites[v]
                    idx = index(alphabet, seq[j])
                    if idx > 0 && !varmask[v, idx] && variants[sites[v]] in keys(AMBIGS)
                        seq_[j] = AMBIGS[variants[sites[v]]]
                    end
                    v += 1
                end
            end
            writeentry(fw, name, join(seq_, ""))
            progress("filtering unsupported variants: processed $i sequences .. ")
        end
    end
    done()
end

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

function main(args)
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

    rates = rategrid(parsed_args["grid-density"], parsed_args["target-rate"])

    const alphabet = "ACGT"
        
    counts = countmsa(parsed_args["msa"], alphabet)

    conditionals, scalers = gridscores(counts, rates, alphabet)

    lls, ws = mcmc(
        conditionals, scalers,
        chain_length=parsed_args["chain-length"], burnin_fraction=parsed_args["burnin-fraction"])

    priors, posteriors = postproc(conditionals', ws)

    variants = callvariants(
        counts, rates, priors, posteriors, alphabet,
        target_rate=parsed_args["target-rate"], posterior_threshold=parsed_args["posterior-threshold"], json_file = parsed_args["json-report"])
        
    if "filter" in keys(parsed_args)
        filterseqs(parsed_args["msa"], variants, parsed_args["filter"], alphabet)
    end
end

main(ARGS)