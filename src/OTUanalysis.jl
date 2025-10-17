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
using DecisionTree

include("constants.jl")
include("load-data.jl")
include("check-data.jl")
include("restructure-data.jl")
include("plot.jl")
include("plsr.jl")
include("postprocessing.jl")
include("random-forest.jl")

# define parameters
env_var = "ST" # The environmental variable to be used for the analysis. Can be "AT", "ST", or "SM".
otu_ids = ["OTU$(lpad(i, 4, '0'))" for i in 1:1] # The OTU IDs to be used for the analysis.
span = 30 * 24 # The time span (in hours) before the sampling date which will be analysed.
season = "all" # The meteorological season which should be included. Can also be "all" to select all seasons.
cdna = false # If true, the data is cDNA.

step = 0.1 # The step size for the selectivity ratio calculation.
n_folds = 50000 # Number of folds for cross-validation.
sig_niveau = 0.05 # The significance level for the selectivity ratio.
smooth = 0.15 # The smoothing parameter for the loess smoothing.
countRange = false # If true, the count range is calculated.

date_col = "datetime" # The name of the column containing the datetime.
id_col = "OTU_ID" # The name of the column containing the OTU IDs.
sampling_date_west = "23.07.2023 11:00" # The date at which the samples were taken in the western study region.
sampling_date_east = "22.07.2023 11:00" # The date at which the samples were taken in the eastern study region.
start_date = "" # e.g."01.01.2022 12:00" # Oprional! The start date of the analysis. If set, span will be ignored. Format: "dd.mm.yyyy HH:MM"
end_date = "" # e.g. "01.01.2023 12:00" # Optional! The end date of the analysis. If not set, the sampling date will be used. Format: "dd.mm.yyyy HH:MM"

plot = true # If true, plots are generated.
plot_pdf = true # If true, plots are saved as pdf.
plot_png = true # If true, plots are saved as png.

saveFrequencies = true # If true, the frequencies (input for PLSR model) are saved.
plot_type = "all" # Can be "all", "points_raw", "points_smoothed", "points_smoothed_with_sig", "line", or "line_sig".
table_form = "vertical" # The form of the output table (.csv). Can be "vertical" or "horizontal".

random_forest = true # If true, a random forest model is used to calculate feature importance.
group_by = "month" # The aggregation for the environmental variables in the random forest model. Can be "month", "year", or "all".

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

# load otu data
if cdna
    df_dna = read_csv("./data/19032025_cDNA_1_clr_sorted.csv")
else
    df_dna = read_csv("./data/19032025_DNA_1_clr_sorted.csv")
end

# plsr
results = Dict()
vec_rmse = []
for otu_id in otu_ids
    results[otu_id] = get_selectivity_ratio!(df_env, df_dna, cdna, otu_id, span, step, n_folds, date_col, id_col, sampling_date_west, sampling_date_east, start_date, end_date, env_var, vec_rmse, season, smooth, plot, plot_pdf, plot_png, countRange, saveFrequencies, plot_type, sig_niveau)
    if random_forest
        random_forest_importance(df_env, df_dna, cdna, otu_id, span, date_col, id_col, sampling_date_west, sampling_date_east, start_date, end_date, season, plot, plot_pdf, plot_png, group_by)
    end
end

# save rmse
df_rmse = DataFrame(
    rmse=vec_rmse,
    id=otu_ids
)
df_rmse.env = fill(env_var, nrow(df_rmse))
df_rmse.season = fill(season, nrow(df_rmse))

const results_path = "./rmse_results.csv"

function add_rmse!(df_rmse::DataFrame)
    if isfile(results_path)
        df_old = CSV.read(results_path, DataFrame)
        df_all = vcat(df_old, df_rmse)
    else
        df_all = df_rmse
    end

    CSV.write(results_path, df_all)
end

add_rmse!(df_rmse)

# postprocessing
process_results(results, env_var, span, season, table_form)

end # module OTUanalysis
