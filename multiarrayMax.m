function indexes = multiarrayMax(A)

    sizeA = size(A); % number of QOI, number of time points, number of parameters
    
    indexes = zeros(1,25);
    maximums = zeros(1,35);

    % for each quantity of interest, determine which parameter is the most influential

    for i = 1:sizeA(1) 
        sens_max_vector = max(A(1,:,:)); % sens_max_vector is 1 by 35. looking at QOI 1, returns the max sensitivity magnitude over the time points for each parameter
        [val idx] = max(sens_max_vector); % returns the highest magnitude of sensitivity seen for any parameter, idx is the index of the param with highest sensitivity
        indexes(i) = idx;
    end

end