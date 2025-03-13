function plot_selectivity_ratio(df::DataFrame, cdna::Bool, otu_id::String, span::Number, env_var::String, season::String="all", save_pdf::Bool=false, save_png::Bool=false, plot_type::String="all")
    """
    Trunctuates a Dataframe with a datetime column to a specific sub-dataframe.

    Parameters:
    - df::DataFrame: Processed DataFrame containing plsr results
    - cdna::Bool: Is the data cDNA? 
    - otu_id::String: ID for the OTUs.
    - span::Number: span (in hours) before the sampling date to which the DataFrame should be trunctuated.  
    - env_var::String: environmental variable to process (either "AT", "ST", "SM")
    - season::String: Optional, the meteorolocical season which should be included. Includes all seasons if not set.
    - save_pdf::Bool: Otional, specifies if the plot should be safed as pdf
    - save_png::Bool: Otional, specifies if the plot should be safed as png
    - plot_type::String: Optional, specifies the type of plot will be generated. Can be "all", "points_raw", "points_smoothed", "points_smoothed_with_sig", "line", or "line_sig"

    Returns:
    - DataFrame: DataFrame containing selectivity ratios (sel_ratio) with significance (significance), p-value (pval), environmental value (x), and smoothed selectivity ratio for plotting (sel_ratio_smooth).
    """
    df_clean = filter(x -> (ismissing(x.sel_ratio) || !isnan(x.sel_ratio)), df)

    if cdna
        cdna_indicator = "c"
    else
        cdna_indicator = ""
    end

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

    title = "$otu_id$(cdna_indicator), for $span hours $season_title"
    if plot_type == "all"
        p = Gadfly.plot(
            layer(
                df_clean,
                x=:x,
                y=:explained_var_smooth_sig,
                color=[colorant"black"],
                Geom.line,
                Theme(line_width=2pt)
            ),
            layer(
                df_clean,
                x=:x,
                y=:explained_var_smooth,
                color=[colorant"grey"],
                Geom.line,
                Theme(line_width=2pt)
            ),
            layer(
                df_clean,
                x=:x,
                y=:explained_var,
                color=:significance,
                Geom.bar
            ),
            Coord.Cartesian(ymin=-1.2, ymax=1.2),
            Scale.color_discrete_manual(palette_sign..., levels=[false, true]),
            Guide.xlabel(x_label),
            Guide.ylabel("Selectivity ratio (smoothed and scaled)"),
            Guide.title(title),
            Theme(background_color="white",
                highlight_width=0mm)
        )
    end

    if plot_type == "points_raw"
        p = Gadfly.plot(
            layer(
                df_clean,
                x=:x,
                y=:explained_var,
                color=:significance,
                Geom.point
            ),
            Coord.Cartesian(ymin=-1.2, ymax=1.2),
            Scale.color_discrete_manual(palette_sign..., levels=[false, true]),
            Guide.xlabel(x_label),
            Guide.ylabel("Selectivity ratio"),
            Guide.title(title),
            Theme(background_color="white",
                highlight_width=0mm)
        )
    end

    if plot_type == "points_smoothed"
        p = Gadfly.plot(
            layer(
                df_clean,
                x=:x,
                y=:explained_var_smooth,
                Geom.point
            ),
            Coord.Cartesian(ymin=-1.2, ymax=1.2),
            Guide.xlabel(x_label),
            Guide.ylabel("Selectivity ratio (smoothed and scaled)"),
            Guide.title(title),
            Theme(
                default_color=palette_sign[1],
                background_color="white",
                highlight_width=0mm)
        )
    end

    if plot_type == "points_smoothed_with_sig"
        p = Gadfly.plot(
            layer(
                df_clean,
                x=:x,
                y=:explained_var_smooth,
                color=:significance,
                Geom.point
            ),
            Coord.Cartesian(ymin=-1.2, ymax=1.2),
            Scale.color_discrete_manual(palette_sign..., levels=[false, true]),
            Guide.xlabel(x_label),
            Guide.ylabel("Selectivity ratio (smoothed and scaled)"),
            Guide.title(title),
            Theme(background_color="white",
                highlight_width=0mm)
        )
    end

    if plot_type == "line"
        p = Gadfly.plot(
            layer(
                df_clean,
                x=:x,
                y=:explained_var_smooth,
                Geom.line,
                Theme(
                    default_color=palette_sign[1],
                    line_width=2pt)
            ),
            Coord.Cartesian(ymin=-1.2, ymax=1.2),
            Guide.xlabel(x_label),
            Guide.ylabel("Selectivity ratio (smoothed and scaled)"),
            Guide.title(title),
            Theme(
                background_color="white",
                highlight_width=0mm)
        )
    end

    if plot_type == "line_sig"
        p = Gadfly.plot(
            layer(
                df_clean,
                x=:x,
                y=:explained_var_smooth_sig,
                Geom.line,
                Theme(
                    default_color=palette_sign[2],
                    line_width=2pt)
            ),
            Coord.Cartesian(ymin=-1.2, ymax=1.2),
            Guide.xlabel(x_label),
            Guide.ylabel("Selectivity ratio (smoothed and scaled)"),
            Guide.title(title),
            Theme(
                default_color=palette_sign[2],
                background_color="white",
                highlight_width=0mm)
        )
    end

    if save_pdf
        draw(PDF("./$(env_var)_$(otu_id)$(cdna_indicator)_$(span)_$(season).pdf", 14cm, 14cm), p)
    end

    if save_png
        draw(PNG("./$(env_var)_$(otu_id)$(cdna_indicator)_$(span)_$(season).png", 14cm, 14cm), p)
    end

    return p
end
