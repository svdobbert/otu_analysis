function check_environmental_input(df::DataFrame, col_name::String, first_date::String, last_date::String)
    """
    Checks if the Dataframes containing the environmental input data have the correct format.

    Parameters:
    - df::DataFrame: Dataframe to be checked.
    - col_name::String: Name of the date column, which should be in the Dataframe
    - first_date::String: First Date which should be included in the Dataframe
    - last_date::String: Last Date which should be included in the Dataframe

    Returns:
    - info if the Dataframe has the expected Format
    """

    if col_name in names(df)
        println("Column $col_name is present in the DataFrame.")
    else
        println("Column $col_name is missing in the DataFrame.")
    end

    first_value = first(df[:, col_name])
    if first_value == first_date
        println("The first value in column $col_name is $first_date")
    else
        println("The first value in column $col_name is not $first_date")
    end

    last_value = last(df[:, col_name])
    if last_value == last_date
        println("The last value in column $col_name is $last_date")
    else
        println("The last value in column $col_name is not $last_date")
    end

    non_float64_nullable_columns = [names(df)[col] for col in 2:ncol(df) if eltype(df[:, col]) != Union{Float64,Missing}]

    if isempty(non_float64_nullable_columns)
        println("All data columns are of type Float64?.")
    else
        println("The following columns are NOT of type Float64?: ", non_float64_nullable_columns)
    end
end

