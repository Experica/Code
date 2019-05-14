function apply_scale(scale::ContinuousColorScale,
                     aess::Vector{Gadfly.Aesthetics}, datas::Gadfly.Data...)
    cmin = Inf
    cmax = -Inf
    for data in datas
        if data.color === nothing
            continue
        end

        for c in data.color
            if c === NA
                continue
            end

            c = convert(Float64, c)
            if c < cmin
                cmin = c
            end

            if c > cmax
                cmax = c
            end
        end
    end

    if cmin == Inf || cmax == -Inf
        return nothing
    end

    if scale.minvalue != nothing
        cmin = scale.minvalue
    end

    if scale.maxvalue  != nothing
        cmax = scale.maxvalue
    end

    cmin, cmax = promote(cmin, cmax)

    ticks, viewmin, viewmax = Gadfly.optimize_ticks(cmin, cmax)
    if ticks[1] == 0 && cmin >= 1
        ticks[1] = 1
    end

    ticks=ticks[cmin.<=ticks.<=cmax]
    #cmin = ticks[1]
    #cmax = ticks[end]
    cspan = cmax != cmin ? cmax - cmin : 1.0

    for (aes, data) in zip(aess, datas)
        if data.color === nothing
            continue
        end

        aes.color = DataArray(RGB{Float32}, length(data.color))
        apply_scale_typed!(aes.color, data.color, scale, cmin, cspan)

        color_key_colors = Dict{Color, Float64}()
        color_key_labels = Dict{Color, String}()

        tick_labels = identity_formatter(ticks)
        for (i, j, label) in zip(ticks, ticks[2:end], tick_labels)
            r = (i - cmin) / cspan
            c = scale.f(r)
            color_key_colors[c] = r
            color_key_labels[c] = label
        end
        c = scale.f((ticks[end] - cmin) / cspan)
        color_key_colors[c] = (ticks[end] - cmin) / cspan
        color_key_labels[c] = tick_labels[end]

        function labeler(xs)
            [get(color_key_labels, x, "") for x in xs]
        end

        aes.color_function = scale.f
        aes.color_label = labeler
        aes.color_key_colors = color_key_colors
        aes.color_key_continuous = true
    end
end
