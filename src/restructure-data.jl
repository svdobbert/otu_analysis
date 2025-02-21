function trunctuate_df(df::DataFrame, date_col::String, sampling_date::String, span::Number)
    """
    Trunctuates a Dataframe with a datetime column to a specific sub-dataframe.

    Parameters:
    - df::DataFrame: DataFrame to be trunctuated.
    - datecol::String: Name of the column containing the datetime.
    - sampling_date::String: Date at which the samples were taken.
    - span::Number: span (in hours) before the sampling date to which the DataFrame should be trunctuated.           

    Returns:
    - DataFrame: The processed DataFrame.
    """
    date_format = "dd.mm.yyyy HH:MM"
    end_date = DateTime(sampling_date, date_format)
    start_date = end_date - Hour(span)

    @info "Including values beween $start_date and $end_date."

    filtered_df = filter(row -> start_date ≤ row[Symbol(date_col)] ≤ end_date, df)

    return filtered_df
end

function is_in_season(date::DateTime, season::String)
    """
    Calculates the meteorolocical season from dates and checks if it is in season.

    Parameters:
    - date::DateTime: A vector of dates.
    - season::String: The season to check for.

    Returns:
    - Vector: a vector of seasons (as strings)
    """

    month = Dates.month(Date.(date))
    return (season == "winter" && month in (12, 1, 2)) ||
           (season == "spring" && month in (3, 4, 5)) ||
           (season == "summer" && month in (6, 7, 8)) ||
           (season == "autumn" && month in (9, 10, 11)) ||
           (season == "all" && month in (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))
end

function count_frequencies(df::DataFrame, date_col::String, range_values)
    """
    Counts frequencies for values above specific values in a range

    Parameters:
    - df::DataFrame: DataFrame with original values.
    - datecol::String: Name of the column containing the datetime.
    - range_values: Vector containing values above which values are counted.

    Returns:
    - DataFrame: The processed DataFrame, containing values as columns and frequencies for each site as rows.
    """

    df_counts = DataFrame()
    range_pairs = [(range_values[i-1], range_values[i+1]) for i in 2:length(range_values)-1]

    for col in names(df)
        if col != date_col
            counts = [count(x -> !ismissing(x) && lower ≤ x < upper, df[!, col])
                      for (lower, upper) in range_pairs]
            df_counts = vcat(df_counts, DataFrame(counts', :auto))
        end
    end

    range_names = [Symbol("$(range_values[i])") for i in 2:length(range_values)-1]

    rename!(df_counts, range_names)

    df_counts.position = names(select(df, Not(Symbol(date_col))))

    return df_counts
end

function prepare_data(df_env::DataFrame, df_otu::DataFrame, otu_id::String, span::Number, step::Number, date_col::String, id_col::String, sampling_date_west::String, sampling_date_east::String, env_var::String, season::String="all")
    """
    Trunctuates a Dataframe with a datetime column to a specific sub-dataframe.

    Parameters:
    - df_env::DataFrame: DataFrame containing environmental data
    - df_otu::DataFrame: DataFrame containing otu data
    - otu_id::String: ID for the OTUs.
    - span::Number: span (in hours) before the sampling date to which the DataFrame should be trunctuated.  
    - step:: Number: step by which the environmental data should be counted.   
    - date_col::String: Name of the column containing the datetime.
    - id_col::String: Name of the column containing the otu ids.
    - sampling_date_west::String: Date at which the samples were taken in the west.
    - sampling_date_east::String: Date at which the samples were taken in the west.
    - env_var::String: environmental variable to process (either "AT", "ST", "SM")
    - season::String: The meteorolocical season which should be included. Can also be "all" to select all seasons.

    Returns:
    - DataFrame: The processed DataFrame.
    """
    if !(id_col in names(df_otu))
        error("The selected id column name ($id_col) is not contained in the dataframe.")
    end

    if !any(occursin(env_var, s) for s in df_env[:, Symbol("env_var")])
        error("The selected env_var ($env_var) is not contained in the dataframe. Available values are $(unique(df_env[:, Symbol("env_var")]))")
    end

    df = filter(row -> row[Symbol("env_var")] == env_var, df_env)
    df = select(df, Not(Symbol("env_var")))

    min_val = minimum(map(x -> minimum(skipmissing(x)), eachcol(select(df, Not(Symbol(date_col))))))
    max_val = maximum(map(x -> maximum(skipmissing(x)), eachcol(select(df, Not(Symbol(date_col))))))
    range_values = collect(min_val:step:max_val)

    @info "Environmental data ranging from $min_val to $max_val."

    date_format = "dd.mm.yyyy HH:MM"

    west_columns = filter(c -> occursin("W", string(c)), names(df))
    east_columns = filter(c -> occursin("E", string(c)), names(df))

    df_west = df[:, west_columns]
    df_west[!, Symbol(date_col)] = DateTime.(df[!, Symbol(date_col)], date_format)

    df_east = df[:, east_columns]
    df_east[!, Symbol(date_col)] = DateTime.(df[!, Symbol(date_col)], date_format)

    df_west_trunctuated = trunctuate_df(df_west, date_col, sampling_date_west, span)
    df_east_trunctuated = trunctuate_df(df_east, date_col, sampling_date_east, span)

    df_west_trunctuated = filter(row -> is_in_season(row[Symbol(date_col)], season), df_west_trunctuated)
    df_east_trunctuated = filter(row -> is_in_season(row[Symbol(date_col)], season), df_east_trunctuated)

    if nrow(df_west_trunctuated) < 1 || nrow(df_east_trunctuated) < 1
        throw(ErrorException("There are now hours within the selected span and season: $season"))
    end

    df_east_frequencies = count_frequencies(df_east_trunctuated, date_col, range_values)
    df_west_frequencies = count_frequencies(df_west_trunctuated, date_col, range_values)

    df_frequencies = vcat(df_west_frequencies, df_east_frequencies)

    selected_row = df_otu[df_otu[!, Symbol(id_col)].==otu_id, :]

    if nrow(selected_row) > 1
        throw(ErrorException("There is more than one row with the given ID: $otu_id"))
    end

    if nrow(selected_row) < 1
        error("Incorrect OTU ID. The selected OTU ID ($otu_id) is not contained in the id column.")
    end

    row_vector = collect(selected_row[1, :])

    cleaned_vector = filter(x -> x isa Number, row_vector)

    if nrow(df_frequencies) != length(cleaned_vector)
        throw(ErrorException("There are no the correct number of numeric values in the otu dataframe"))
    end

    df_otu_selected = DataFrame()
    df_otu_selected.values = cleaned_vector
    df_otu_selected.position = names(select(df_otu, Not(Symbol(id_col))))

    df_plsr_input = innerjoin(df_frequencies, df_otu_selected, on=:position)

    return df_plsr_input
end
