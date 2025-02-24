function plot_selectivity_ratio(df::DataFrame, otu_id::String, span::Number, env_var::String, season::String="all", save_pdf::Bool=false, save_png::Bool=false)
    """
    Trunctuates a Dataframe with a datetime column to a specific sub-dataframe.

    Parameters:
    - df::DataFrame: Processed DataFrame containing plsr results
    - otu_id::String: ID for the OTUs.
    - span::Number: span (in hours) before the sampling date to which the DataFrame should be trunctuated.  
    - env_var::String: environmental variable to process (either "AT", "ST", "SM")
    - season::String: Optional, the meteorolocical season which should be included. Includes all seasons if not set.
    - save_pdf::Bool: Otional, specifies if the plot should be safed as pdf
    - save_png::Bool: Otional, specifies if the plot should be safed as png

    Returns:
    - DataFrame: DataFrame containing selectivity ratios (sel_ratio) with significance (significance), p-value (pval), environmental value (x), and smoothed selectivity ratio for plotting (sel_ratio_smooth).
    """
    df_clean = filter(x -> (ismissing(x.sel_ratio_smooth) || !isnan(x.sel_ratio_smooth)), df)

    if !any(occursin(env_var, s) for s in ["AT", "ST", "SM"])
        @warn "The selected env_var ($env_var) is invalid. Available values are $(["AT", "ST", "SM"])"
        x_label = "Environmental variable"
    end

    if env_var == "AT"
        x_label = "Air Temperature [°C]"
    end

    if env_var == "ST"
        x_label = "Soil Temperature [°C]"
    end

    if env_var == "SM"
        x_label = "Soil Moisture [m²/m²]"
    end

    if season == "all"
        season_title = ""
    else
        season_title = ", season = $season"
    end

    title = "$otu_id, for $span hours $season_title"

    p = Gadfly.plot(
        layer(
            df_clean,
            x=:x,
            y=:sel_ratio_smooth,
            color=:significance,
            Geom.point
        ),
        Scale.color_discrete_manual(palette_sign...),
        Guide.xlabel(x_label),
        Guide.ylabel("Selectivity Ratio (smoothed)"),
        Guide.title(title),
        Theme(background_color="white",
            highlight_width=0mm)
    )

    if save_pdf
        draw(PDF("./$(env_var)_$(otu_id)_$span_$(season).pdf", 14cm, 14cm), p)
    end

    if save_png
        draw(PNG("./$(env_var)_$(otu_id)_$(span)_$(season).png", 14cm, 14cm), p)
    end

    return p
end
