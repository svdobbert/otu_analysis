function transform_data(df::DataFrame, env_var::String, span::Number, date_col::String, sampling_date_west::String, sampling_date_east::String, start_date::String, end_date::String, group_by::String, season::String="all")
    """
    Prepares the DataFrame, trunctuating and transforming the data

    Parameters:
    - df::DataFrame: DataFrame containing environmental data
    - env_var::String: environmental variable to process (either "AT", "ST", "SM")
    - span::Number: span (in hours) before the sampling date to which the DataFrame should be trunctuated.  
    - date_col::String: Name of the column containing the datetime.
    - sampling_date_west::String: Date at which the samples were taken in the west.
    - sampling_date_east::String: Date at which the samples were taken in the east.
    - start_date::String: Optional, the start date of the analysis. If set, span will be ignored.
    - end_date::String: Optional, the end date of the analysis. If not set, the sampling date will be used.
    - group_by::String: The grouping of the data. Can be "month", "year", or "all".
    - season::String: The meteorolocical season which should be included. Can also be "all" to select all seasons.

    Returns:
    - DataFrame: The processed DataFrame.
    """
    df = filter(row -> row[Symbol("env_var")] == env_var, df)
    df = select(df, Not(Symbol("env_var")))
    df_trunctuated = prepare_env_data(df, span, date_col, sampling_date_west, sampling_date_east, start_date, end_date, season)

    if group_by == "month"
        df_trunctuated.group = Dates.month.(df_trunctuated[!, Symbol(date_col)])
        df_grouped = combine(groupby(df_trunctuated, :group), names(df_trunctuated, Not(Symbol(date_col), :group)) .=> mean)
        group = unique(string.(sort(df_trunctuated.group)))

        month_map = Dict(
            "1" => "Jan", "2" => "Feb", "3" => "Mar",
            "4" => "Apr", "5" => "May", "6" => "Jun",
            "7" => "Jul", "8" => "Aug", "9" => "Sep",
            "10" => "Oc", "11" => "Nov", "12" => "Dec"
        )
        group = [month_map[m] for m in group]
    end

    if group_by == "year"
        df_trunctuated.group = Dates.year.(df_trunctuated[!, Symbol(date_col)])
        df_grouped = combine(groupby(df_trunctuated, :group), names(df_trunctuated, Not(Symbol(date_col), :group)) .=> mean)
        group = unique(string.(df_trunctuated.group))
    end

    if group_by == "all"
        df_trunctuated.group = fill("all", nrow(df_trunctuated))
        df_grouped = combine(groupby(df_trunctuated, :group), names(df_trunctuated, Not(Symbol(date_col), :group)) .=> mean)
        group = unique(string.(df_trunctuated.group))
    end
    df_transposed = DataFrame(permutedims(df_grouped))[2:end, :]

    rename!(df_transposed, group)

    return df_transposed
end

function random_forest_importance(df_env::DataFrame, df_otu::DataFrame, cdna::Bool, otu_id::String, span::Number, date_col::String, id_col::String, sampling_date_west::String, sampling_date_east::String, start_date::String, end_date::String, season::String="all", plot=true, plot_pdf::Bool=false, plot_png::Bool=false, group_by::String="all")
    """
    Performs a random forest regression to determine feature importance of environmental variables for a given OTU.
        
    Parameters:
    - df_env::DataFrame: DataFrame containing environmental data
    - df_otu::DataFrame: DataFrame containing otu data
    - cdna::Bool: Is the data cDNA? 
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
    selected_row = df_normalized[df_normalized[!, Symbol(id_col)].==otu_id, :]

    if nrow(selected_row) > 1
        throw(ErrorException("There is more than one row with the given ID: $otu_id"))
    end

    if nrow(selected_row) < 1
        error("Incorrect OTU ID. The selected OTU ID ($otu_id) is not contained in the id column.")
    end

    row_vector = collect(selected_row[1, :])

    cleaned_vector = filter(x -> x isa Number, row_vector)

    if nrow(df_transformed) != length(cleaned_vector)
        throw(ErrorException("There are not the correct number of numeric values in the otu dataframe"))
    end

    df_otu_selected = DataFrame()
    df_otu_selected.values = cleaned_vector
    df_otu_selected.position = names(select(df_otu, Not(Symbol(id_col))))

    df_rf_input = innerjoin(df_transformed, df_otu_selected, on=:position)
    df_rf_input = dropmissing(select(df_rf_input, Not(:position)))

    X = df_rf_input[:, Not(:values)]
    x_values = names(X)
    X = Matrix{Float64}(X)
    y = convert(Vector{Float64}, df_rf_input.values)

    #  normalizing X and y 
    std_X = std(X, dims=1)
    std_X[std_X.==0] .= 1
    X = (X .- mean(X, dims=1)) ./ std_X

    std_y = std(y, dims=1)
    std_y = std_y == 0 ? 1 : std_y
    y = (y .- mean(y, dims=1)) ./ std_y

    # train random forest classifier
    model = RandomForestRegressor(n_trees=100, max_depth=10, n_subfeatures=3)
    # fit model and get feature importances
    model = DecisionTree.fit!(model, X, y)
    importances = impurity_importance(model)
    feature_names = x_values
    df_importance = DataFrame(feature=feature_names, importance=importances)
    sort!(df_importance, :importance, rev=true)

    # if sum(df_importance.importance) != 1
    #     throw(ErrorException("The sum of the importances is not equal to 1. Please check the model."))
    # end

    if cdna
        cdna_indicator = "c"
    else
        cdna_indicator = ""
    end

    CSV.write("./$(otu_id)$(cdna_indicator)_$(span)_$(season)_feature_importance.csv", df_importance)

    # plot feature importance
    if plot
        if season == "all"
            season_title = ""
        else
            season_title = ", season = $season"
        end

        title = "$otu_id$(cdna_indicator), for $span hours $season_title"
        p = Gadfly.plot(
            df_importance,
            x=:feature,
            y=:importance,
            Geom.bar,
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
            draw(PDF("./$(otu_id)$(cdna_indicator)_$(span)_$(season)_feature_importance.pdf", 14cm, 14cm), p)
        end

        if plot_png
            draw(PNG("./$(otu_id)$(cdna_indicator)_$(span)_$(season)_feature_importance.png", 14cm, 14cm), p)
        end

        display(p)
        return p
    end
end
