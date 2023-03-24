
### log pdf of a multivariate normal
function logpfun_mv(x, μ, σ, w)
    logw = log.(w)
    K = length(logw)
    N = size(x, 2)
    sum(logsumexp(logw[k] + logpdf(MvNormal(μ[k], σ), x[:,n]) for k in 1:K) for n in 1:N)
end


@model function nuc_pair_model(x, nc, nsites, σ,  p,  wc=ones(nc)/nc, wp=ones(nsites-1)/(nsites-1), w=wc.*wp', q = 2)
    
    n = size(x, 2)
    startpos ~ filldist(Uniform(-1500, 1500), nc)
    delta    ~ filldist(Gamma(q*p, 1/q), nc, nsites)
    configs = startpos .+ cumsum(delta, dims=2)
    
    μ = [[configs[i, j], configs[i, j+1]] for i = 1:nc, j=1:(nsites-1)]
    
    #x ~ filldist(MixtureModel(Normal.(vec(μ), w)), n)
    Turing.@addlogprob! logpfun_mv(x, vec(μ), σ, w)
end



@model function smf_bp_single(reads)

    open_p   ~ Beta(1, 1)
    closed_p ~ Beta(1, 1)
    
    footprint_p ~ filldist(Beta(1, 1), size(reads, 1))
    footprint ~ Bernoulli.(footprint_p)

    for i = 1:size(reads, 1)
        reads[i, :] ~ filldist(Bernoulli(footprint[i] ? open_p : closed_p), size(reads, 2))
    end
end


@model function smf_bp_single_nv(reads)

    open_p   ~ Beta(1, 1)
    closed_p ~ Beta(1, 1)
    
    footprint_p ~ filldist(Beta(1, 1), size(reads, 1))

    footprint = Vector{Bool}(undef, size(reads, 1))
    for i = 1:size(reads, 1)
        footprint[i] ~ Bernoulli(footprint_p[i])
    end
    
    for j = 1:size(reads, 2), i = 1:size(reads, 1)
        reads[i, j] ~ Bernoulli(footprint[i] ? open_p : closed_p)
    end

end

function smf_bp_total(total_meth, n)
    open_p   ~ Beta(1, 1)
    closed_p ~ Beta(1, 1)
    m = length(total_meth)

    footprint_p ~ filldist(Beta(1, 1), m)
    footprint ~ Bernoulli.(footprint_p)

    total_meth ~ filldist(Binomial(n, footprint ? open_p : closed_p), m)
    
end
