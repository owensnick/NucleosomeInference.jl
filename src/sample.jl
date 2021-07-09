

function randclasses(classconfig, n, σ=10, q=2, seed=1618)
    Random.seed!(seed)

    nc = length(classconfig)

    configs = Vector{Float64}[]

    for i = 1:nc
        delta = rand(Gamma(q*classconfig[i].p, 1/q), classconfig[i].m)
        config = classconfig[i].start .+ cumsum(delta)
        push!(configs, config)
    end
    
    ln = length.(configs) .- 1
    class = rand(1:nc, n)    
    left = [rand(1:ln[c]) for c in class]
    right = left .+ 1
    left_obs = [rand(Normal(configs[c][l], σ)) for (c, l) in zip(class, left)]
    right_obs = [rand(Normal(configs[c][l], σ)) for (c, l) in zip(class, right)]
    
    x = [left_obs' ; right_obs']
    
    obs   =  DataFrame(class=class, left=left, right=right, left_obs=left_obs, right_obs=right_obs)
    

    
    (configs=configs, n=n, obs=obs, classconfig=classconfig, x=x)

end

function simpletwoclass(n=1000, σ=10)
    classconfig = [(start=-900, p=180, m=6), (start=50, p=180, m=6)]
    randclasses(classconfig, n, σ)
end

