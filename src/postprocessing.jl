function process_results(results::Dict, span::Number, step::Number, season::String, sort::String=vertical)
    """
    Restructures the results of the previous analysis and saves the to a csv file

    Parameters:
    - results::Dict: The results of the previous analysis.
    - span::Number: span (in hours) before the sampling date to which the DataFrame should be trunctuated.  
    - step::Number: step by which the environmental data should be counted.   
    - season::String: The meteorolocical season which should be included. Can also be "all" to select all seasons.
    - sort::String: The sorting of the data. Can be either "vertical" or "horizontal".

    Returns:
    - DataFrame: The processed DataFrame.
    """

    if sort === "horizontal"
        all_x_values = unique(vcat([df.x for df in values(results)]...))
        sort!(all_x_values)

        df_full = DataFrame(x=all_x_values)

        # Perform outer join iteratively
        df_merged = df_full
        for (otu_id, result) in results
            df_result = DataFrame()
            df_result[!, Symbol("$(otu_id)_sel_ratio")] = result.sel_ratio
            df_result[!, Symbol("$(otu_id)_p_val")] = result.p_val
            df_result[!, Symbol("$(otu_id)_significance")] = result.significance
            df_result[!, Symbol("$(otu_id)_sel_ratio_smooth")] = result.sel_ratio_smooth
            df_result[!, Symbol("$(otu_id)_explained_var")] = result.explained_var
            df_result[!, Symbol("$(otu_id)_explained_var_smooth")] = result.explained_var_smooth
            df_result[!, :x] = result.x
            df_merged = outerjoin(df_merged, df_result, on=:x)
        end

        df_merged .= coalesce.(df_merged, NaN)
    end

    if sort === "vertical"
        df_merged = DataFrame()
        for (otu_id, result) in results
            df_result = DataFrame(result)
            df_result.OTU_ID = fill(otu_id, nrow(df_result))
            df_merged = vcat(df_merged, df_result)
        end

        df_merged .= coalesce.(df_merged, NaN)
    end

    formatted_time = Dates.format(now(), "yyyy-mm-dd_HH:MM")
    CSV.write("./$(formatted_time)_$(env_var)_$(span)_$(season).csv", df_merged)
end