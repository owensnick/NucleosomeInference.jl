
function sample_mh_jump(model, samples=1000 ; δ=180, σ=5, p=0.2)
    mh = MH(:startpos => x -> [MixtureModel(Normal.([v -  δ, v + δ], σ), [p/2, 1-p, p/2]) for v in x])
    sample(model, mh, samples)
end