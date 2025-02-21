function read_csv(file_path::String; output_path::Union{String,Nothing}=nothing)
    """
    Reads a CSV file where numbers use commas as decimal separators and converts them to proper Float64 format.

    Parameters:
    - file_path::String: Path to the input CSV file.
    - output_path::Union{String, Nothing}: Optional path to save the modified file. 
                                           If not provided, it will process in-memory.

    Returns:
    - DataFrame: The processed DataFrame with numbers properly formatted.
    """
    file_content = read(file_path, String)

    file_content = replace(file_content, r"(\d),(\d)" => s"\1.\2")

    if output_path !== nothing
        open(output_path, "w") do f
            write(f, file_content)
        end
        df = CSV.read(output_path, DataFrame)
        return filter(row -> any(!ismissing, row), df)
    else
        io = IOBuffer(file_content)
        df = CSV.read(io, DataFrame)
        return filter(row -> any(!ismissing, row), df)
    end
end

