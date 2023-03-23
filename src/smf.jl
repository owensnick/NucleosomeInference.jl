"""
    sample_smf_data(tf_footprints, seq, nreads=10, open_meth_prob=0.9, closed_meth_prob=0.05, fwprob=0.5)

    Sample methylation data from a sequence with a given number of reads, and a given probability of methylation for open and closed sites.

"""
function sample_smf_data(tf_footprints, seq, nreads=10, open_meth_prob=0.9, closed_meth_prob=0.05, fwprob=0.5)

    pos_A = seq .== DNA_A
    pos_T = seq .== DNA_T

    open_A = seq .== DNA_A
    open_T = seq .== DNA_T

    for (s, e) in tf_footprints
        open_A[s:e] .= false
        open_T[s:e] .= false
    end
    

    fA = findall(pos_A)
    fT = findall(pos_T)

    sampled_reads = Vector{Int}[]

    for i = 1:nreads
        fr = rand(Bernoulli(fwprob))

        if fr
            r = sample_pos(fA, open_A, open_meth_prob, closed_meth_prob)
        else
            r = sample_pos(fT, open_T, open_meth_prob, closed_meth_prob)
        end
        push!(sampled_reads, r)
    end

    sampled_reads
end

"""
    sample_pos(pos, open, open_meth_prob, closed_meth_prob)

    Sample methylation sites at a given set of positions, with a given probability of methylation for open and closed sites.
"""
function sample_pos(pos, open, open_meth_prob, closed_meth_prob)
    positions = Int[]

    for p in pos
        if open[p]
            if rand(Bernoulli(open_meth_prob))
                push!(positions, p)
            end
        else
            if rand(Bernoulli(closed_meth_prob))
                push!(positions, p)
            end
        end
    end

    positions
end