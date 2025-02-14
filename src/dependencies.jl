# dependencies

using Pkg

Pkg.activate(".")

Pkg.add("Gadfly")               # a system for plotting and visualization 
Pkg.add("Arrow")                # Arrow storage and file format
Pkg.add("CategoricalArrays")    # similar to the factor type in R
Pkg.add("CSVFiles")             # read/write CSV and similar formats
Pkg.add("Downloads")            # file downloads
Pkg.add("DataFrames")           # versatile tabular data format
Pkg.add("GZip")                 # utilities for compressed files
Pkg.add("Tar")                  # tar archive utilities
Pkg.add("MultivariateStats")    # for multivariate statistics and data analysis
Pkg.add("Plots")                # powerful convenience for visualization in Julia
Pkg.add("MLJ")                  # a Machine Learning Framework for Julia
Pkg.add("Dates")                # provides types for working with dates
Pkg.add("Statistics")           # basic statistics functionality.
Pkg.add("Missings")             # Convenience functions for working with missing values in Julia
Pkg.add("LinearAlgebra")        # common and useful linear algebra operations 
Pkg.add("Compose")              # a declarative vector graphics system written in Julia
Pkg.add("ColorSchemes")         # a set of pre-defined ColorSchemes
Pkg.add("Jchemo")               # Chemometrics and machine learning on high-dimensional data
