function plot_genes_cell_segments(df_spatial::DataFrame, polygons::Array{Matrix{Float64}, 1}=Matrix{Float64}[]; 
        genes::Union{Vector, Symbol, String}="Podxl", 
        colors::Union{Vector, Symbol, String}="blue", 
        bg_color::Union{Vector, Symbol, String}="gray90",
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
