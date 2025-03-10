module OTUanalysis

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

include("constants.jl")
include("load-data.jl")
include("check-data.jl")
include("restructure-data.jl")
include("plot.jl")
include("plsr.jl")

# define parameters
env_var = "ST"
otu_ids = ["OTU1000", "OTU1001"]
span = 365 * 24
step = 0.1
date_col = "datetime"
id_col = "OTU_ID"
sampling_date_west = "23.07.2023 11:00"
sampling_date_east = "22.07.2023 11:00"
season = "all"
smooth = 0.5
plot = true
plot_pdf = false
plot_png = false
countRange = false
saveFrequencies = true

# load environmental data
df_at = read_csv("./data/AT15_Data_2009_2023_fixed.csv")
df_st = read_csv("./data/ST15_Data_2009_2023.csv")
df_sm = read_csv("./data/SM15_2009_2023.csv")

df_at.env_var = fill("AT", nrow(df_at))
df_st.env_var = fill("ST", nrow(df_st))
df_sm.env_var = fill("SM", nrow(df_sm))
df_env = vcat(df_at, df_st, df_sm)

last_date_west = "23.07.2023 11:00"
last_date_east = "22.07.2023 11:00"

# check environmental data
check_environmental_input(df_at, "datetime", "15.09.2009 01:00", "23.07.2023 11:00")
check_environmental_input(df_st, "datetime", "15.09.2009 01:00", "23.07.2023 11:00")
check_environmental_input(df_sm, "datetime", "15.09.2009 01:00", "23.07.2023 11:00")

# load otu data
df_dna = read_csv("./data/clr_sorted_DNA_OTU_PLSR_final.csv")
df_vdna = read_csv("./data/clr_sorted_cDNA_OTU_PLSR_final.csv")
cdna = true

results = Dict()
for otu_id in otu_ids
    results[otu_id] = get_selectivity_ratio(df_env, df_dna, cdna, otu_id, span, step, date_col, id_col, sampling_date_west, sampling_date_east, env_var, season, smooth, plot, plot_pdf, plot_png, countRange, saveFrequencies)
end

end # module OTUanalysis
