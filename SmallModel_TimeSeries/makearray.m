function parameter_array = makearray(parameter_structure)

    parameter_cell = struct2cell(parameter_structure);
    parameter_array = cell2mat(parameter_cell);

end