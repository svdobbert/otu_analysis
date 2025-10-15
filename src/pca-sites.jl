module PCAsites

using XLSX
using DataFrames
using PlotlyJS          # interface to the [plotly.js] visualization library
using Tar               # tar archive utilities
using MLJ
PCA_model = @load PCA pkg = "MultivariateStats" verbosity=0
using MultivariateStats
using CSV


function drop_cols_with_missing(df::DataFrame, threshold::Float64)
    keep = [name for name in names(df) if count(ismissing, df[!, name]) / nrow(df) ≤ threshold]
    return df[:, keep]
end

type = "DNA"
time = "01d"

function pca(type::String)
    """
    Calculates a PCA and plots the results using PlotlyJS.

    Parameters:
    - `type::String`: The type of data to be used for the PCA. Can be "cDNA" or "DNA".

    Returns:
    - A 3d PCA-Plot.
    """
    filename = "data/$(type)_mean_sorted_for_PCA.xlsx"
    filename_taxa = "data/$(type)_taxa.csv"
    filename_colors = "colours_PCA_order.xlsx"

    df = DataFrame(XLSX.readtable(filename, "Sheet1", infer_eltypes=true))
    # df = drop_cols_with_missing(df, 0.8)

    # df_clean = dropmissing(df)
    df_clean = coalesce.(df, 0)
    features = Symbol.(string.(df_clean.x))
    df_clean = df_clean[:, Not([1, ncol(df_clean)])]
    df_clean = DataFrame([[names(df_clean)]; collect.(eachrow(df_clean))], [:column; Symbol.(axes(df_clean, 1))])
    rename!(df_clean, names(df_clean)[1] => :OTU_ID)

    df_taxa = CSV.read(filename_taxa, DataFrame)
    df_taxa = leftjoin(df_clean, df_taxa, on = :OTU_ID)
    df_clean = hcat(df_taxa[:, Symbol("Order")], df_clean)

    for c in names(df_clean[!, 3:end])
        if eltype(df_clean[!, c]) <: AbstractString
            df_clean[!, c] = parse.(Float64, df_clean[!, c])
        end
    end
    model = PCA_model(maxoutdim=3)
    mach = machine(model, df_clean[!, 3:end])
    fit!(mach)

    components = MLJ.transform(mach, df_clean[!, 3:end])
    components.otu_id = df_clean[!, 2]
    components.order = df_clean[!, 1]
    total_var = report(mach).tprincipalvar / report(mach).tvar
    projection = fitted_params(mach).projection
    loadings = projection' .* report(mach).principalvars
	df_loadings = DataFrame(loadings, :auto)
	
	rename!(df_loadings, features)

	CSV.write("./pca/pca_loadings_$(type)_sites.csv", df_loadings)

    # plot PCA
    df_colors = DataFrame(XLSX.readtable(filename_colors, "Tabelle1", infer_eltypes=true))
    df_colors = rightjoin(df_colors, components, on = :order)[!, [:colour]]
    colors_points = df_colors.colour
    colors_points = coalesce.(colors_points, "#cdcdbf") 
    components.color = colors_points

    groups = unique(components.order)
    scatter_points = [
        scatter3d(
            x = components.x1[components.order .== g],
            y = components.x2[components.order .== g],
            z = components.x3[components.order .== g],
            mode = "markers",
            marker = attr(
                size = 3,
                color = components.color[components.order .== g],
                ),
            name = string(g),
            text = components.otu_id[components.order .== g],
            hoverinfo = "text"
        )
        for g in groups
    ]

    arrow_traces = [
        scatter3d(
            x=[0, loadings[1, i]],
            y=[0, loadings[2, i]],
            z=[0, loadings[3, i]],
            mode="lines+text",
            line=attr(color="#646464", width=2),
            text=[nothing, features[i]],
            textposition="top center",
            textfont=attr(
                size=12, 
                color="#646464"), 
            name=features[i],
            showlegend=false
        )
        for i in 1:length(features)
    ]

    p = PlotlyJS.plot(
        [scatter_points; arrow_traces...],
        Layout(
            title="$(type) sites, total variance: $(round(total_var, digits=3))",
            scene=attr(
                xaxis=attr(title="PC1"),
                yaxis=attr(title="PC2"),
                zaxis=attr(title="PC3")
            )
        )
    )

    filename = "pca_$(type)_sites.html"
    savefig(p, "./pca/$(filename)")
end

pca(type)

end # module PCAsites
