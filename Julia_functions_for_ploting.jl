# This is a modified function based on Baysor's plot_molecules that allows visualizing the expression of the selected genes on top of the cell segmentation.
function plot_genes_cell_segments(df_spatial::DataFrame, polygons::Array{Matrix{Float64}, 1}=Matrix{Float64}[]; 
        genes::Union{Vector, Symbol, String}="Podxl", 
        colors::Union{Vector, Symbol, String}="blue", 
        bg_color::Union{Vector, Symbol, String}="gray95",
        genesize =2, segline_size=0.5, 
        canvas_size=(5000,6000),x_lims=nothing, y_lims=nothing, 
        ann::Union{<:AbstractVector, Symbol, Nothing} = nothing,
        ann_colors::Union{Nothing, Dict} = nothing, 
        show_legend=(ann !== nothing), legend_fontsize=12, transparency=0.5
    )
    other_genes=unique(df_spatial.gene_id[Not(in.(df_spatial.gene_id, [Set(genes)]))])
    other_colors=repeat([bg_color],length(other_genes))
    all_genes=[genes; other_genes]
    all_colors=[colors; other_colors]
    map_color=Dict(all_genes .=> all_colors)
    df_spatial=transform(df_spatial, :gene_id => ByRow(x -> map_color[x]) => :new_color)
    fig=B.plot_molecules(df_spatial, polygons, 
        color=:new_color, markersize=genesize, 
        size=canvas_size,poly_strokewidth=segline_size, 
        xlims=x_lims, ylims=y_lims,
        annotation=ann, ann_colors=ann_colors,
        legend=show_legend, fontsize=legend_fontsize, alpha=transparency
    )
    return fig
end

function plot_genes_overlaid_cells(df_spatial::DataFrame, polygons::Array{Matrix{Float64}, 1}=Matrix{Float64}[]; 
        genes::Union{Vector, Symbol, String}="Podxl", 
        colors::Union{Vector, Symbol, String}="blue", 
        bg_color::Union{Vector, Symbol, String}="gray95",
        markersize =2, segline_size=0.5, offset=(0, 0),
        canvas_size=(5000,6000),x_lims=nothing, y_lims=nothing, 
        annotation::Union{<:AbstractVector, Symbol, Nothing}=nothing,
        ann_colors::Union{Nothing, Dict}=nothing, noise_ann = nothing,
        show_legend=false,legend_fontsize=12, transparency=0.5,
        legend_kwargs::Union{Dict, NamedTuple, Nothing}=nothing
    )
    legend_args_default = (bgcolor=Colors.RGBA(1, 1, 1, 0.85),);
    legend_kwargs = B.update_args(legend_args_default, legend_kwargs)
    if annotation !== nothing
        if typeof(annotation) === Symbol
            annotation = df_spatial[!,annotation]
        end
        annotation = ["$a" for a in annotation]
        if noise_ann !== nothing
            noise_ann = "$noise_ann"
        end
    end
    other_genes=unique(df_spatial.gene_id[Not(in.(df_spatial.gene_id, [Set(genes)]))])
    other_colors=repeat([(bg_color,0.1)],length(other_genes))
    all_genes=[genes; other_genes]
    all_colors=[colors; other_colors]
    map_color=Dict(all_genes .=> all_colors)
    df_spatial=transform(df_spatial, :gene_id => ByRow(x -> map_color[x]) => :new_color)
    xlims = something(x_lims, B.val_range(df_spatial.x))
    ylims = something(y_lims, B.val_range(df_spatial.y))
    fig = MK.Figure(resolution=canvas_size)
    fig[1, 1] = MK.Axis(fig; xticklabelsize=12, yticklabelsize=12, xticksvisible=false, xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false);
    ann_vals = annotation[annotation .!= noise_ann] |> unique |> sort
    c_map = Colors.distinguishable_colors(length(ann_vals), Colors.colorant"#007a10", lchoices=range(20, stop=70, length=15))
    for (color, ann) in zip(c_map, ann_vals)
        style_dict = (ann_colors === nothing) ? Dict() : Dict(:color => ann_colors[ann])
        MK.scatter!(df_spatial.x[annotation .== ann] .+ offset[1], df_spatial.y[annotation .== ann] .+ offset[2];
            strokewidth=0, markersize=markersize, label=ann, color=color, style_dict...)
    end
    if show_legend
        MK.axislegend(;legend_kwargs...)
    end
    MK.poly!([MK.Point2.(eachrow(p .+ [offset[1] offset[2]])) for p in polygons]; strokecolor="black", color="transparent", strokewidth=segline_size, label="")
    colors2 = df_spatial[!,:new_color]
    MK.scatter!(df_spatial.x .+ offset[1], df_spatial.y .+ offset[2]; color=colors2,
            strokewidth=0, markersize=markersize)
    MK.xlims!(MK.current_axis(), xlims .+ offset[1])
    MK.ylims!(MK.current_axis(), ylims .+ offset[2])
    return MK.current_figure()
end
