### plot sampled smf data
function plot_smf_sample(data; c=[:red, :black], kwargs...)

    # plot the data
    p = plot(; kwargs...)


    plot_reads!(data.fA, data.forward_reads, 1; c=c[1])
    plot_reads!(data.fT, data.reverse_reads, size(data.forward_reads, 2) + 1; c=c[2])

    ### mark footprints with vertical line
    for fp in data.footprints
        vline!([first(fp), last(fp)], c=:red, ls=:dash, lab="")
    end
    p

end


### plot set of reads
function plot_reads!(pos, reads, s=1; c=:black)

    # plot the data
    p = plot!()
    
    for i ∈ axes(reads, 2)
        x = pos[reads[:, i]]
        pi = i + s - 1
        plot!([minimum(x), maximum(x)], [pi, pi], lab="", c=c)
        scatter!(x, fill(pi, length(x)), lab="", marker=(:vline, 5), c=c)
    end

    p
end

