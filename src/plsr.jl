function remove_constant_columns(df::DataFrame)
    non_constant_columns = [length(unique(skipmissing(df[:, col]))) > 1 for col in names(df)]
    return df[:, non_constant_columns]
end

function get_selectivity_ratio(df::DataFrame, df_otu::DataFrame, otu_id::String, span::Number, step::Number, date_col::String, id_col::String, sampling_date_west::String, sampling_date_east::String, env_var::String, season::String="all")
    """
    Trunctuates a Dataframe with a datetime column to a specific sub-dataframe.

    Parameters:
    - df::DataFrame: DataFrame containing environmental data
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
    - DataFrame: DataFrame containing selectivity ratios (sel_ratio) with significance (sign), p-value (pval), and environmental value (x).
    """

    df_processed = prepare_data(df, df_otu, otu_id, span, step, date_col, id_col, sampling_date_west, sampling_date_east, env_var, season)
    df_cleaned = dropmissing(select(df_processed, Not(:position)))

    X = remove_constant_columns(df_cleaned[:, Not(:values)])
    env_values = names(X)
    X = Matrix{Float64}(X)
    y = convert(Vector{Float64}, df_cleaned.values)

    if size(X, 2) == 0
        error("All features were removed due to being constant. Check your data preprocessing.")
    end

    #  normalizing X and y 
    std_X = std(X, dims=1)
    std_X[std_X.==0] .= 1  # Replace zero standard deviations with 1 to prevent division by zero
    X = (X .- mean(X, dims=1)) ./ std_X

    std_y = std(y, dims=1)
    std_y = std_y == 0 ? 1 : std_y  # Prevent division by zero for y
    y = (y .- mean(y, dims=1)) ./ std_y

    # Ensure no NaNs or Infs
    X[isnan.(X).|isinf.(X)] .= 0
    y[isnan.(y).|isinf.(y)] .= 0

    # Cross Validation parameters
    n_folds = 10
    n_samples = size(X, 1)
    n_permutations = 1000  # Number of permutations
    n_features = size(X, 2)
    subset_size = round(Int, 0.8 * n_samples)
    folds = [rand(1:n_samples, subset_size) for _ in 1:n_folds]

    # Step 1: Calculate actual selectivity ratios
    selectivity_ratios_all_folds = []
    coefficient_signs_all_folds = []
    p_values = zeros(n_features)

    n_features_original = length(env_values)

    for fold_idx in 1:n_folds
        test_indices = folds[fold_idx]
        train_indices = setdiff(1:n_samples, test_indices)

        # Split data into train and test sets
        X_train = X[train_indices, :]
        y_train = y[train_indices]

        # Perform feature selection (remove near-zero variance features)
        variances = var(X_train, dims=1)
        keep_features = vec(variances .> 1e-6)
        X_train_selected = X_train[:, keep_features]

        println("Remaining features after selection for fold $fold_idx: ", sum(keep_features))

        global nlv_value = min(3, size(X_train, 1) - 1)  
        if nlv_value < 1
            error("Not enough samples for PLS (nlv = $nlv_value).")
        end

        # Fit the PLS model (ensure it works with reduced features)
        pls_model = Jchemo.plswold(X_train_selected, y_train, nlv=nlv_value, scal=true)

        Wx = pls_model.W
        feature_variance = var(X_train_selected, dims=1)

        # Avoid division by zero
        feature_variance[feature_variance.==0] .= 1

        explained_variance = sum(Wx .^ 2, dims=2)
        selectivity_ratios = explained_variance ./ feature_variance'

        sign_of_loadings = sign.(Wx)
        coefficient_signs = sign.(mean(sign_of_loadings, dims=2))

        # Create a vector of selectivity ratios for all features (matching original feature set)
        selectivity_ratios_full = fill(NaN, n_features_original)
        coefficient_signs_full = fill(NaN, n_features_original)

        # Assign selectivity ratios to the selected features
        selectivity_ratios_full[keep_features.==1] .= selectivity_ratios
        coefficient_signs_full[keep_features.==1] .= coefficient_signs

        # Store the selectivity ratios for this fold
        push!(coefficient_signs_all_folds, coefficient_signs_full)
        push!(selectivity_ratios_all_folds, selectivity_ratios_full)
    end

    # Concatenate the selectivity ratios across folds (they are now aligned with the same features)
    selectivity_ratios_matrix = hcat(selectivity_ratios_all_folds...)
    println("Number of NaN values in the selectivity ratios matrix: ", count(isnan, selectivity_ratios_matrix))
    selectivity_ratios_df = DataFrame(selectivity_ratios_matrix, :auto)
    selectivity_ratios_mean = [mean(skipmissing(row)) for row in eachrow(selectivity_ratios_df)]

    coefficient_signs_matrix = hcat(coefficient_signs_all_folds...)
    coefficient_signs_mean = sign.(mean(coefficient_signs_matrix, dims=2))

    # Step 2: Permutation test for p-values
    p_values = zeros(n_features)

    for feature_idx in 1:n_features
        permuted_ratios = Float64[]

        for _ in 1:n_permutations
            permuted_y = shuffle(y)

            pls_model_perm = Jchemo.plswold(X, permuted_y, nlv=nlv_value, scal=true)

            Wx_perm = pls_model_perm.W
            explained_variance_perm = sum(Wx_perm .^ 2, dims=2)
            feature_variance_perm = var(X, dims=1)

            permuted_selectivity_ratios = explained_variance_perm ./ feature_variance_perm'

            push!(permuted_ratios, mean(permuted_selectivity_ratios))  # Take mean across folds
        end
        # Calculate p-value: Proportion of permutations greater than or equal to the actual value
        p_values[feature_idx] = mean(abs.(permuted_ratios) .>= abs(selectivity_ratios_mean[feature_idx]))
    end

    selectivity_ratios_with_sign = vec(coefficient_signs_mean .* selectivity_ratios_mean)

    plsr_result = DataFrame(
        sel_ratio=vec(selectivity_ratios_with_sign),
        p_val=p_values,
        sign=p_values .< 0.1,
        x=env_values
    )

    return plsr_result
end
