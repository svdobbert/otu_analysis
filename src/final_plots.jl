using Arrow
using CategoricalArrays
using CSVFiles
using Downloads
using DataFrames
using Plots
using Dates
using Gadfly
using Compose
using ColorSchemes
using StatsBase
using CSV
using Cairo
using Dates

include("load-data.jl")
include("constants.jl")

# plot rmse
df_rmse = read_csv("./rmse_results.csv")
sort!(df_rmse, :rmse, rev=false)
df_rmse.x = df_rmse.id .* "_" .* df_rmse.env .* "_" .* df_rmse.season

p = Gadfly.plot(
    df_rmse,
    x=:x,
    y=:rmse,
    Geom.bar,
    Guide.xlabel("Model Input"),
    Guide.ylabel("RMSE"),
    Coord.cartesian(ymin=0, ymax=1),
    Theme(
        default_color=main_col,
        background_color="white",
        highlight_width=0mm)
)

formatted_time = Dates.format(now(), "yyyy-mm-dd_HH_MM")
draw(PDF("./$(formatted_time)_rmse.pdf", 14cm, 14cm), p)
draw(PNG("./$(formatted_time)_rmse.png", 14cm, 14cm), p)

display(p)
