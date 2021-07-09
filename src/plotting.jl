

function plotconfig(config, dx=10 ; kwargs...)
    ex = extrema(config.x)
    bins = first(ex):dx:last(ex)
    stephist(vec(config.x), bins=bins, group=repeat(config.obs.class, inner=2), leg=:outertopright; kwargs...)
end