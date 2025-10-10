module RandomForest

using Arrow
using CategoricalArrays
using CSVFiles
using Downloads
using DataFrames
using GZip
using Tar
using MultivariateStats
using Plots
using MLJ
using Dates
using Statistics
using Missings
using LinearAlgebra
using Gadfly
using Compose
using ColorSchemes
using Jchemo
using DataFrames
using StatsBase
using CSV
using Loess
using Cairo
using Dates
using DecisionTree


include("constants.jl")
include("load-data.jl")
include("check-data.jl")
include("restructure-data.jl")
include("random-forest-plsr.jl")

otu_ids = ["OTU$(lpad(i, 4, '0'))" for i in 1:2] # The OTU IDs to be used for the analysis.
span = 30 * 24 # The time span (in hours) before the sampling date which will be analysed.
season = "all" # The meteorological season which should be included. Can also be "all" to select all seasons.
date_col = "datetime" # The name of the column containing the datetime.
id_col = "OTU_ID" # The name of the column containing the OTU IDs.
sampling_date_west = "23.07.2023 11:00" # The date at which the samples were taken in the western study region.
sampling_date_east = "22.07.2023 11:00" # The date at which the samples were taken in the eastern study region.
start_date = "" # e.g."01.01.2022 12:00" # Oprional! The start date of the analysis. If set, span will be ignored. Format: "dd.mm.yyyy HH:MM"
end_date = "" # e.g. "01.01.2023 12:00" # Optional! The end date of the analysis. If not set, the sampling date will be used. Format: "dd.mm.yyyy HH:MM"
plot_png = true # If true, plots are saved as png.
group_by = "year" # The aggregation for the environmental variables in the random forest model. Can be "month", "year", or "all".

# load environmental data
df_at = read_csv("./data/AT15_Data_2009_2023_fixed.csv")
df_st = read_csv("./data/ST15_Data_2009_2023.csv")
df_sm = read_csv("./data/SM15_2009_2023.csv")

df_at.env_var = fill("AT", nrow(df_at))
df_st.env_var = fill("ST", nrow(df_st))
df_sm.env_var = fill("SM", nrow(df_sm))
df_env = vcat(df_at, df_st, df_sm)

# check environmental data
check_environmental_input(df_at, "datetime", "15.09.2009 01:00", "23.07.2023 11:00")
check_environmental_input(df_st, "datetime", "15.09.2009 01:00", "23.07.2023 11:00")
check_environmental_input(df_sm, "datetime", "15.09.2009 01:00", "23.07.2023 11:00")

df_cdna = read_csv("./data/19032025_cDNA_1_clr_sorted.csv")
df_dna = read_csv("./data/19032025_DNA_1_clr_sorted.csv")

df_cdna.type = fill("cDNA", nrow(df_cdna))
df_dna.type = fill("DNA", nrow(df_dna))
df_otu = vcat(df_dna, df_cdna)

function random_forest_importance_complete(df_env::DataFrame, df_otu::DataFrame, otu_id::String, span::Number, date_col::String, id_col::String, sampling_date_west::String, sampling_date_east::String, start_date::String, end_date::String, season::String="all", group_by::String="all")
    """
    Performs a random forest regression to determine feature importance of environmental variables for a given OTU.
        
    Parameters:
    - df_env::DataFrame: DataFrame containing environmental data
    - df_otu::DataFrame: DataFrame containing otu data
    - otu_id::String: ID for the OTUs.
    - span::Number: span (in hours) before the sampling date to which the DataFrame should be trunctuated.  
    - date_col::String: Name of the column containing the datetime.
    - id_col::String: Name of the column containing the otu ids.
    - sampling_date_west::String: Date at which the samples were taken in the west.
    - sampling_date_east::String: Date at which the samples were taken in the east.
    - start_date::String: Optional, the start date of the analysis. If set, span will be ignored.
    - end_date::String: Optional, the end date of the analysis. If not set, the sampling date will be used.
    - season::String: Optional, the meteorolocical season which should be included. Includes all seasons if not set.
    - plot::Bool: Otional, specifies if a plot should be returned.
    - plot_pdf::Bool: Otional, specifies if the plot should be safed as pdf
    - plot_png::Bool: Otional, specifies if the plot should be safed as png
    - group_by::String: The grouping of the data. Can be "month", "year", or "all".

    Returns:
    - Feature importance as .csv and (optional) plot.
    """
    df_at_transformed = transform_data(df_env, "AT", span, date_col, sampling_date_west, sampling_date_east, start_date, end_date, group_by, season)
    rename!(df_at_transformed, Dict(names(df_at_transformed)[i] => "AT_" * string(names(df_at_transformed)[i]) for i in 1:length(names(df_at_transformed))))
    df_st_transformed = transform_data(df_env, "ST", span, date_col, sampling_date_west, sampling_date_east, start_date, end_date, group_by, season)
    rename!(df_st_transformed, Dict(names(df_st_transformed)[i] => "ST_" * string(names(df_st_transformed)[i]) for i in 1:length(names(df_st_transformed))))
    df_sm_transformed = transform_data(df_env, "SM", span, date_col, sampling_date_west, sampling_date_east, start_date, end_date, group_by, season)
    rename!(df_sm_transformed, Dict(names(df_sm_transformed)[i] => "SM_" * string(names(df_sm_transformed)[i]) for i in 1:length(names(df_sm_transformed))))

    df_transformed = hcat(df_at_transformed, df_st_transformed, df_sm_transformed)
    df_transformed.position = names(select(df_env, Not(Symbol(date_col), Symbol("env_var"))))

    df_normalized = normalize_otu_data(df_otu, id_col)
    selected_row_dna = df_normalized[(df_normalized[!, Symbol(id_col)].==otu_id).&(df_normalized[!, :type].=="DNA"), :]
    selected_row_cdna = df_normalized[(df_normalized[!, Symbol(id_col)].==otu_id).&(df_normalized[!, :type].=="cDNA"), :]

    if nrow(selected_row_dna) > 1
        throw(ErrorException("There is no DNA data for the given ID: $otu_id"))
    end

    if nrow(selected_row_cdna) > 1
        throw(ErrorException("There is no cDNA data for the given ID: $otu_id"))
    end

    if nrow(selected_row_dna) < 1 | nrow(selected_row_cdna) < 1
        error("Incorrect OTU ID. The selected OTU ID ($otu_id) is not contained in the id column.")
    end

    row_vector_dna = collect(selected_row_dna[1, :])
    row_vector_cdna = collect(selected_row_cdna[1, :])

    cleaned_vector_dna = filter(x -> x isa Number, row_vector_dna)
    cleaned_vector_cdna = filter(x -> x isa Number, row_vector_cdna)

    if nrow(df_transformed) != length(cleaned_vector_dna) | nrow(df_transformed) != length(cleaned_vector_cdna)
        throw(ErrorException("There are not the correct number of numeric values in the otu dataframe"))
    end

    df_otu_selected = DataFrame()
    df_otu_selected.dna_values = cleaned_vector_dna
    df_otu_selected.cdna_values = cleaned_vector_cdna
    df_otu_selected.position = names(select(df_otu, Not([Symbol(id_col), Symbol("type")])))

    df_rf_input = innerjoin(df_transformed, df_otu_selected, on=:position)
    df_rf_input = dropmissing(select(df_rf_input, Not(:position)))

    X = df_rf_input[:, Not([:dna_values, :cdna_values])]
    x_values = names(X)
    X = Matrix{Float64}(X)
    y_dna = convert(Vector{Float64}, df_rf_input.dna_values)
    y_cdna = convert(Vector{Float64}, df_rf_input.cdna_values)

    #  normalizing X and y 
    std_X = std(X, dims=1)
    std_X[std_X.==0] .= 1
    X = (X .- mean(X, dims=1)) ./ std_X

    std_y_dna = std(y_dna, dims=1)
    std_y_dna = std_y_dna == 0 ? 1 : std_y_dna
    y_dna = (y_dna .- mean(y_dna, dims=1)) ./ std_y_dna

    std_y_cdna = std(y_cdna, dims=1)
    std_y_cdna = std_y_cdna == 0 ? 1 : std_y_cdna
    y_cdna = (y_cdna .- mean(y_cdna, dims=1)) ./ std_y_cdna

    # train random forest classifier
    model_dna = RandomForestRegressor(n_trees=100, max_depth=10, n_subfeatures=3)
    model_cdna = RandomForestRegressor(n_trees=100, max_depth=10, n_subfeatures=3)
    # fit model and get feature importances
    model_dna = DecisionTree.fit!(model_dna, X, y_dna)
    model_cdna = DecisionTree.fit!(model_cdna, X, y_cdna)
    importances_dna = impurity_importance(model_dna)
    importances_cdna = impurity_importance(model_cdna)
    feature_names = x_values
    df_importance = DataFrame(
        feature=repeat(feature_names, 2),
        importance=[importances_dna; importances_cdna],
        type=vcat(fill("DNA", length(feature_names)), fill("cDNA", length(feature_names))),
        otu_id=fill(otu_id, 2 * length(feature_names)),
        env=repeat(split.(feature_names, "_") .|> x -> x[1], 2),
        variable=repeat(split.(feature_names, "_") .|> x -> join(x[2:end], "_"), 2)
    )
    sort!(df_importance, :importance, rev=true)

    return df_importance
end

function plot_feature_importance(rf_result::DataFrame, otu_id::String, season::String, span::Int64, plot_pdf::Bool, plot_png::Bool)
    """
    Plots the feature importance for a given OTU.

    Parameters:
    - rf_result::DataFrame: DataFrame containing the feature importance results from random_forest_importance_complete
    - otu_id::String: ID for the OTUs.
    - season::String: The meteorolocical season which should be included. Can also be "all" to select all seasons.
    - span::Number: span (in hours) before the sampling date to which the DataFrame should be trunctuated.  
    - plot_pdf::Bool: Otional, specifies if the plot should be safed as pdf
    - plot_png::Bool: Otional, specifies if the plot should be safed as png 
    Returns:
    - Plot of feature importance.
    """
    df = rf_result[rf_result.otu_id.==otu_id, :]

    if season == "all"
        season_title = ""
    else
        season_title = ", season = $season"
    end

    title = "$otu_id - Feature importance (Random Forest) for $span hours $season_title"
    p = Gadfly.plot(
        df,
        x=:feature,
        y=:importance,
        color=:type,
        Geom.bar(position=:dodge),
        Scale.color_discrete_manual(palette_sign..., levels=["DNA", "cDNA"]),
        Guide.xlabel("Feature"),
        Guide.ylabel("Importance"),
        Guide.title(title),
        Coord.cartesian(ymin=0, ymax=1),
        Theme(
            default_color=main_col,
            background_color="white",
            highlight_width=0mm)
    )

    if plot_pdf
        draw(PDF("./$(otu_id)_$(span)_$(season)_feature_importance.pdf", 14cm, 14cm), p)
    end

    if plot_png
        draw(PNG("./$(otu_id)_$(span)_$(season)_feature_importance.png", 14cm, 14cm), p)
    end

    display(p)
    return p

end

dfs = []
for otu_id in otu_ids
    df = random_forest_importance_complete(df_env, df_otu, otu_id, span, date_col, id_col, sampling_date_west, sampling_date_east, start_date, end_date, season, group_by)
    push!(dfs, df)
end
rf_result = vcat(dfs...)
CSV.write("./$(span)_$(season)_feature_importance.csv", rf_result)

for otu_id in otu_ids
    plot_feature_importance(rf_result, otu_id, season, span, true, true)
end
end # module RandomForest