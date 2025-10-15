module PCAresults

using XLSX
using DataFrames
using PlotlyJS          # interface to the [plotly.js] visualization library
using Tar               # tar archive utilities
using CSV
using MLJ
PCA_model = @load PCA pkg = "MultivariateStats" verbosity=0
using MultivariateStats

function drop_cols_with_missing(df::DataFrame, threshold::Float64)
    keep = [name for name in names(df) if count(ismissing, df[!, name]) / nrow(df) ≤ threshold]
    return df[:, keep]
end

type = "DNA"
time = "01d"
reduce_env = true
digits_T = 0
digits_SM = 2

function pca(type::String, time::String, reduce_env::Bool=false, digits_T::Int=0, digits_SM::Int=2)
    """
    Calculates a PCA and plots the results using PlotlyJS.

    Parameters:
    - `type::String`: The type of data to be used for the PCA. Can be "cDNA" or "DNA".
 	- `time::String`: The time point for which the PCA should be calculated. E.g. "01d", "07d", etc.

    Returns:
    - A 3d PCA-Plot.
    """
    filename_at = "data/$(type)results/$(type)_AT_$(time)_postprocessing_R_explVar.xlsx"
    filename_st = "data/$(type)results/$(type)_ST_$(time)_postprocessing_R_explVar.xlsx"
    filename_sm = "data/$(type)results/$(type)_SM_$(time)_postprocessing_R_explVar.xlsx"
    filename_taxa = "data/$(type)_taxa.csv"
    filename_colors = "colours_PCA_order.xlsx"

    df_at = DataFrame(XLSX.readtable(filename_at, "Sheet1", infer_eltypes=true))

    df_st = DataFrame(XLSX.readtable(filename_st, "Sheet1", infer_eltypes=true))

    df_sm = DataFrame(XLSX.readtable(filename_sm, "Sheet1", infer_eltypes=true))

    df_at.env_var = fill("AT", nrow(df_at))
    df_st.env_var = fill("ST", nrow(df_st))
    df_sm.env_var = fill("SM", nrow(df_sm))
	rename!(df_at, names(df_at)[1] => :x)
	rename!(df_st, names(df_st)[1] => :x)
    rename!(df_sm, names(df_sm)[1] => :x)
    df = vcat(df_at, df_st, df_sm)
    # df = drop_cols_with_missing(df, 0.8)

    # df_clean = dropmissing(df)
    df_clean = coalesce.(df, 0)
    features = Symbol.(string.(df_clean.x) .* "_" .* df_clean.env_var)
    rounded_features = ifelse.(
        df_clean.env_var .== "SM",
        round.(df_clean.x, digits=digits_SM),
        round.(df_clean.x, digits=digits_T),
    )
    rounded_features = Symbol.(string.(rounded_features) .* "_" .* df_clean.env_var)
    env_var = Symbol.(df_clean.env_var)
    df_clean = df_clean[:, Not([1, ncol(df_clean)])]
    df_clean = DataFrame([[names(df_clean)]; collect.(eachrow(df_clean))], [:column; Symbol.(axes(df_clean, 1))])

    df_taxa = CSV.read(filename_taxa, DataFrame)
    df_clean = hcat(df_taxa[:, Symbol("Order")], df_clean)

    # load and fit PCA
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
    unique_names = unique(rounded_features)
    env_var_reduced = Symbol.([split(string(name), "_", limit=2)[2] for name in unique_names])

    reduced_loadings = zeros(size(loadings, 1), length(unique_names))

    for (i, name) in enumerate(unique_names)
        indices = findall(==(name), rounded_features)
        reduced_loadings[:, i] = mean(loadings[:, indices], dims=2)
    end

	CSV.write("./pca/pca_loadings_$(type)_$(time)_explained_variance.csv", df_loadings)

    # plot PCA
    unique_groups = unique(env_var)
    palette = ["#b4b4b4", "#646464", "#000000"]
    group_colors = Dict(g => color for (g, color) in zip(unique_groups, palette))

    if reduce_env
        println("PCA with reduced environmental variables")
        colors_env = [group_colors[env_var_reduced[i]] for i in 1:length(unique_names)]
    else
        println("PCA with all environmental variables")
        colors_env = [group_colors[env_var[i]] for i in 1:length(features)]
    end

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

    plot_loadings = reduce_env ? reduced_loadings : loadings
    plot_features = reduce_env ? unique_names : features

    arrow_traces = [
        scatter3d(
            x=[0, plot_loadings[1, i]],
            y=[0, plot_loadings[2, i]],
            z=[0, plot_loadings[3, i]],
            mode="lines+text",
            line=attr(color=colors_env[i], width=4),
            text=[nothing, plot_features[i]],
            textposition="top center",
            textfont=attr(
                size=14, 
                color=colors_env[i]), 
            name=features[i],
            showlegend=false
        )
        for i in 1:length(plot_features)
    ]

    p = PlotlyJS.plot(
        [scatter_points; arrow_traces...],
        Layout(
            title="$(type) $(time) explained variance, total variance: $(round(total_var, digits=3))",
            scene=attr(
                xaxis=attr(title="PC1"),
                yaxis=attr(title="PC2"),
                zaxis=attr(title="PC3")
            )
        )
    )

    filename = "pca_$(type)_$(time)_explained_variance.html"
    savefig(p, "./pca/$(filename)")
end

pca(type, time, reduce_env, digits_T, digits_SM)

end # module PCAresults
