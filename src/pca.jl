using XLSX
using DataFrames
using PlotlyJS          # interface to the [plotly.js] visualization library
using Tar               # tar archive utilities
using MLJ
using CSV

function drop_cols_with_missing(df::DataFrame, threshold::Float64)
    keep = [name for name in names(df) if count(ismissing, df[!, name]) / nrow(df) ≤ threshold]
    return df[:, keep]
end

type = "DNA"
time = "01d"

function pca(type::String, time::String)
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

    df_at = DataFrame(XLSX.readtable(filename_at, "Sheet1", infer_eltypes=true))
    df_at.x = round.(df_at.x, digits=1)
    df_at = combine(groupby(df_at, :x), names(df_at) .=> mean; keepkeys=false)

    df_st = DataFrame(XLSX.readtable(filename_st, "Sheet1", infer_eltypes=true))
    df_st.x = round.(df_st.x, digits=1)
    df_st = combine(groupby(df_st, :x), names(df_st) .=> mean; keepkeys=false)

    df_sm = DataFrame(XLSX.readtable(filename_sm, "Sheet1", infer_eltypes=true))
    df_sm.x = round.(df_sm.x, digits=2)
    df_sm = combine(groupby(df_sm, :x), names(df_sm) .=> mean; keepkeys=false)

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
    features = Symbol.(string.(round.(df_clean.x, digits = 2)) .* "_" .* df_clean.env_var)
    env_var = Symbol.(df_clean.env_var)
    df_clean = df_clean[:, Not([1, ncol(df_clean)])]
    df_clean = DataFrame([[names(df_clean)]; collect.(eachrow(df_clean))], [:column; Symbol.(axes(df_clean, 1))])


    # load and fit PCA
    PCA = @load PCA pkg = "MultivariateStats"
    mach = machine(PCA(maxoutdim=3), df_clean[!, 2:end])
    fit!(mach)

    components = MLJ.transform(mach, df_clean[!, 2:end])
    components.otu_id = df_clean[!, 1]
    total_var = report(mach).tprincipalvar / report(mach).tvar
    projection = fitted_params(mach).projection
    loadings = projection' .* report(mach).principalvars
	df_loadings = DataFrame(loadings, :auto)
	
	rename!(df_loadings, features)

	CSV.write("./pca/pca_loadings_$(type)_$(time)_explained_variance.csv", df_loadings)

    # plot PCA
    unique_groups = unique(env_var)
    palette = ["#954650", "#9A5554", "#506F8A"]
    group_colors = Dict(g => color for (g, color) in zip(unique_groups, palette))

    colors_env = [group_colors[env_var[i]] for i in 1:length(features)]

    scatter_points = scatter3d(
        x=components.x1, y=components.x2, z=components.x3,
        mode="markers", marker=attr(size=2, color="#A89294"),
        name="OTUs", showlegend=false,
        text=components.otu_id,
        hoverinfo="text",
    )

    arrow_traces = [
        scatter3d(
            x=[0, loadings[1, i]],
            y=[0, loadings[2, i]],
            z=[0, loadings[3, i]],
            mode="lines+text",
            line=attr(color=colors_env[i], width=4),
            text=[nothing, features[i]],
            textposition="top center",
            textfont=attr(size=14, color=colors_env[i]), # adjust text size
            name=features[i],
            showlegend=false
        )
        for i in 1:length(features)
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

pca(type, time)


