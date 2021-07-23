
### log pdf of a multivariate normal
function logpfun_mv(x, μ, σ, w)
    logw = log.(w)
    K = length(logw)
    N = size(x, 2)
    sum(logsumexp(logw[k] + logpdf(MvNormal(μ[k], σ), x[:,n]) for k in 1:K) for n in 1:N)
end


@model function (x, nc, nsites, σ,  p,  wc=ones(nc)/nc, wp=ones(nsites-1)/(nsites-1), w=wc.*wp', q = 2)
    
    n = size(x, 2)
    startpos ~ filldist(Uniform(-1500, 1500), nc)
    delta    ~ filldist(Gamma(q*p, 1/q), nc, nsites)
    configs = startpos .+ cumsum(delta, dims=2)
    
    μ = [[configs[i, j], configs[i, j+1]] for i = 1:nc, j=1:(nsites-1)]
    
    #x ~ filldist(MixtureModel(Normal.(vec(μ), w)), n)
    Turing.@addlogprob! logpfun_mv(x, vec(μ), σ, w)
end