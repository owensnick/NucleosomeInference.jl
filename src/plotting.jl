

function plotconfig(config, dx=10 ; kwargs..)
    xrange = extrema(config.x)
    bins = first(xrange):dx:last(xrange)
    stephist(vec(config.x), bins=bins, group=repeat(config.obs.class, innner=2), leg=:outertopright; kwargs...)
end