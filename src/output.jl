function output_results(dict_results, file_path)
    indentation = 2 # whitespace characters used for indentation
    open(file_path,"w") do io
        JSON.print(io, dict_results, indentation)
    end
end
