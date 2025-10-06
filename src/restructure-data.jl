function trunctuate_df(df::DataFrame, date_col::String, sampling_date::String, span::Number, start_date::String, end_date::String)
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
    if end_date == ""
        try
            end_date = DateTime(sampling_date, date_format)
        catch e
            error("Error: The sampling date '$sampling_date' does not match the format '$date_format'")
        end
    else
        try
            end_date = DateTime(end_date, date_format)
        catch e
            error("Error: The proved end date '$end_date' does not match the format '$date_format'")
        end
        if end_date > DateTime(sampling_date, date_format)
            error("Error: The provided end date '$end_date' is after the sampling date '$sampling_date'")
        end
    end

    if start_date == ""
        start_date = end_date - Hour(span)
    else
        try
            start_date = DateTime(start_date, date_format)
        catch e
            println("Error: The proved start date '$start_date' does not match the format '$date_format'")
            return nothing
        end
    end

    @info "Including values beween $start_date and $end_date."

    filtered_df = filter(row -> start_date ≤ row[Symbol(date_col)] ≤ end_date, df)

    return filtered_df
end


function is_in_season(date::DateTime, season::String)
    """
    Calculates the meteorolocical season from dates.

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


function count_frequencies_range(df::DataFrame, date_col::String, range_values)
    """
    Counts frequencies around specific values in a range

    Parameters:
    - df::DataFrame: DataFrame with original values.
    - datecol::String: Name of the column containing the datetime.
    - range_values: Vector containing values around which values are counted.

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


function count_frequencies(df::DataFrame, date_col::String, range_values)
    """
    Counts frequencies above/below specific values in a range

    Parameters:
    - df::DataFrame: DataFrame with original values.
    - datecol::String: Name of the column containing the datetime.
    - range_values: Vector containing values above/below which values are counted.

    Returns:
    - DataFrame: The processed DataFrame, containing values as columns and frequencies for each site as rows.
    """

    df_counts = DataFrame()

    for col in names(df)
        if col != date_col
            counts = [count(x -> !ismissing(x) && (value ≥ 0 ? x ≥ value : x ≤ value), df[!, col])
                      for value in range_values]
            df_counts = vcat(df_counts, DataFrame(counts', :auto))
        end
    end

    range_names = [Symbol("$(range_values[i])") for i in 1:length(range_values)]

    rename!(df_counts, range_names)

    df_counts.position = names(select(df, Not(Symbol(date_col))))

    return df_counts
end


function prepare_env_data(df::DataFrame, span::Number, date_col::String, sampling_date_west::String, sampling_date_east::String, start_date::String, end_date::String, season::String="all")
    """
    Prepares the DataFrame, trunctuating data with different sampling dates seperately to a specific time span before the sampling date.

    Parameters:
    - df::DataFrame: DataFrame containing environmental data
    - span::Number: span (in hours) before the sampling date to which the DataFrame should be trunctuated.  
    - date_col::String: Name of the column containing the datetime.
    - sampling_date_west::String: Date at which the samples were taken in the west.
    - sampling_date_east::String: Date at which the samples were taken in the west.
    - start_date::String: Start date of the analysis. If set, span will be ignored.
    - end_date::String: End date of the analysis. If not set, the sampling date will be used.
    - season::String: The meteorolocical season which should be included. Can also be "all" to select all seasons.

    Returns:
    - DataFrame: The processed DataFrame.
    """
    date_format = "dd.mm.yyyy HH:MM"

    west_columns = filter(c -> occursin("W", string(c)), names(df))
    east_columns = filter(c -> occursin("E", string(c)), names(df))

    df_west = df[:, west_columns]
    df_west[!, Symbol(date_col)] = DateTime.(df[!, Symbol(date_col)], date_format)

    df_east = df[:, east_columns]
    df_east[!, Symbol(date_col)] = DateTime.(df[!, Symbol(date_col)], date_format)

    df_west_trunctuated = trunctuate_df(df_west, date_col, sampling_date_west, span, start_date, end_date)
    df_east_trunctuated = trunctuate_df(df_east, date_col, sampling_date_east, span, start_date, end_date)

    df_west_trunctuated = filter(row -> is_in_season(row[Symbol(date_col)], season), df_west_trunctuated)
    df_east_trunctuated = filter(row -> is_in_season(row[Symbol(date_col)], season), df_east_trunctuated)

    if nrow(df_west_trunctuated) < 1 || nrow(df_east_trunctuated) < 1
        throw(ErrorException("There are no hours within the selected span and season: $season"))
    end

    df_east_trunctuated = select(df_east_trunctuated, Not(Symbol(date_col)))
    df_trunctuated = hcat(df_west_trunctuated, df_east_trunctuated)

    return df_trunctuated
end

function normalize_otu_data(df::DataFrame, id_col::String)
    """
    Normalizes and cleans a DataFrame, keeping non-numeric columns unmodified.

    Parameters:
    - df::DataFrame: DataFrame containing otu data
    - id_col::String: Name of the column containing the otu ids.
    
    Returns:
    - DataFrame: The processed DataFrame.
    """
    numeric_cols = [c for c in names(df) if c != id_col && eltype(df[!, c]) <: Number]
    non_numeric_cols = [c for c in names(df) if c != id_col && !(eltype(df[!, c]) <: Number)]

    Y = Matrix{Float64}(select(df, numeric_cols))
    std_Y = std(Y, dims=1)
    Y = (Y .- mean(Y, dims=1)) ./ std_Y
    df_Y = DataFrame(Y, Symbol.(numeric_cols))

    for c in non_numeric_cols
        df_Y[!, Symbol(c)] = df[!, c]
    end
    df_Y[!, Symbol(id_col)] = df[!, Symbol(id_col)]

    return df_Y
end

function prepare_data(df_env::DataFrame, df_otu::DataFrame, cdna::Bool, otu_id::String, span::Number, step::Number, date_col::String, id_col::String, sampling_date_west::String, sampling_date_east::String, start_date::String, end_date::String, env_var::String, season::String="all", countRange::Bool=true, saveFrequencies::Bool=true)
    """
    Trunctuates a Dataframe with a datetime column to a specific sub-dataframe and counts frequencies for a following plsr analysis.

    Parameters:
    - df_env::DataFrame: DataFrame containing environmental data
    - df_otu::DataFrame: DataFrame containing otu data
    - cdna::Bool: Is the data cDNA? 
    - otu_id::String: ID for the OTUs.
    - span::Number: span (in hours) before the sampling date to which the DataFrame should be trunctuated.  
    - step:: Number: step by which the environmental data should be counted.   
    - date_col::String: Name of the column containing the datetime.
    - id_col::String: Name of the column containing the otu ids.
    - sampling_date_west::String: Date at which the samples were taken in the west.
    - sampling_date_east::String: Date at which the samples were taken in the west.
    - start_date::String: Start date of the analysis. If set, span will be ignored.
    - end_date::String: End date of the analysis. If not set, the sampling date will be used.
    - env_var::String: environmental variable to process (either "AT", "ST", "SM")
    - season::String: The meteorolocical season which should be included. Can also be "all" to select all seasons.
    - countRange::Bool: Should the frequencies be counted as range around the values or above/below the values.
    - saveFrequencies::Bool: Should the frequencies table be safed.

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

    df_trunctuated = prepare_env_data(df, span, date_col, sampling_date_west, sampling_date_east, start_date, end_date, season)

    west_columns = filter(c -> occursin("W", string(c)), names(df))
    east_columns = filter(c -> occursin("E", string(c)), names(df))

    date_format = "dd.mm.yyyy HH:MM"

    df_east_trunctuated = df_trunctuated[:, east_columns]
    df_east_trunctuated[!, Symbol(date_col)] = df_trunctuated[!, Symbol(date_col)]
    df_west_trunctuated = df_trunctuated[:, west_columns]
    df_west_trunctuated[!, Symbol(date_col)] = df_trunctuated[!, Symbol(date_col)]

    if countRange
        df_east_frequencies = count_frequencies_range(df_east_trunctuated, date_col, range_values)
        df_west_frequencies = count_frequencies_range(df_west_trunctuated, date_col, range_values)
    else
        df_east_frequencies = count_frequencies(df_east_trunctuated, date_col, range_values)
        df_west_frequencies = count_frequencies(df_west_trunctuated, date_col, range_values)
    end

    df_frequencies = vcat(df_west_frequencies, df_east_frequencies)

    # normalize OTU-data
    df_Y = normalize_otu_data(df_otu, id_col)

    selected_row = df_Y[df_Y[!, Symbol(id_col)].==otu_id, :]

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

    if cdna
        cdna_indicator = "c"
    else
        cdna_indicator = ""
    end

    if saveFrequencies
        CSV.write("./$(env_var)_$(otu_id)$(cdna_indicator)_$(span)_$(season)_frequencies.csv", df_plsr_input)
    end

    return df_plsr_input
end
