

### function to load fasta files in the folder data
### records description is formatted: hg19_dna range=chr10:71108610-71108717 5'pad=0 3'pad=0 strand=+ repeatMasking=non

function loadseq(file="hk1_reg_region.fa", projdir=getprojectdir())
    filepath = joinpath(projdir, "data", file)
    reader = open(FASTA.Reader, filepath) 
    record = first(collect(reader))
    close(reader)
    seq = LongDNA{4}(sequence(record))

    desc = split(replace(description(record), "range=" => ""))[2]
    chromloc = split(desc, r"[:-]")
    (seq=seq, chrom=first(chromloc), loc=parse(Int, chromloc[2]):parse(Int, chromloc[3]))

end


function sample_hk1(;nreads=10, open_meth_prob=0.9, closed_meth_prob=0.05, fwprob=0.5)
    seq = loadseq()
    ### look for footprint sequences within seq
    ### TATTCCA for nfat
    ### GTACTCGA for NKX2-2
    ### TGTTTGT for foxa2
    footprint_nfat  = findall(ExactSearchQuery(dna"TATTCCA"), seq.seq)
    footprint_nkx2  = findall(ExactSearchQuery(dna"GTACTCGA"), seq.seq)
    footprint_foxa2 = findall(ExactSearchQuery(dna"TGTTTGT"), seq.seq)
    
    ### check that only one footprint sequence is found for each of nfat, nkx2, foxa2
    @assert length(footprint_nfat) == 1
    @assert length(footprint_nkx2) == 1
    @assert length(footprint_foxa2) == 1
    
    footprints = reduce(vcat, [footprint_nfat, footprint_nkx2, footprint_foxa2])
    footprint_labels = ["NFAT", "NKX2-2", "FOXA2"]
 
    sample = sample_smf_data(footprints, seq.seq, nreads, open_meth_prob, closed_meth_prob, fwprob)

    (sample..., seq=seq, footprints=footprints, labels=footprint_labels, open_meth_prob=open_meth_prob, closed_meth_prob=closed_meth_prob, fwprob=fwprob)    
end

"""
    sample_smf_data(tf_footprints, seq, nreads=10, open_meth_prob=0.9, closed_meth_prob=0.05, fwprob=0.5)

    Sample methylation data from a sequence with a given number of reads, and a given probability of methylation for open and closed sites.

"""
function sample_smf_data(tf_footprints, seq, nreads=10, open_meth_prob=0.9, closed_meth_prob=0.05, fwprob=0.5)

  

    pos_A = seq .== DNA_A
    pos_T = seq .== DNA_T
   
    fA = findall(pos_A)
    fT = findall(pos_T)

    a_in_fp = [!any(fp -> f ∈ fp, tf_footprints) for f in fA]
    t_in_fp = [!any(fp -> f ∈ fp, tf_footprints) for f in fT]

    
    fn = rand(Binomial(nreads, fwprob))
    rn = nreads - fn

    forward_reads = mapreduce(_ -> rand.(Bernoulli.(ifelse.(a_in_fp, open_meth_prob, closed_meth_prob))), hcat, 1:fn)
    reverse_reads = mapreduce(_ -> rand.(Bernoulli.(ifelse.(t_in_fp, open_meth_prob, closed_meth_prob))), hcat, 1:rn)

    
    (fA=fA, fT=fT, a_in_fp=a_in_fp, t_in_fp, forward_reads=forward_reads, reverse_reads=reverse_reads)
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