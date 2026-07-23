function N = getNorms(A)
N = zeros(25,35,7); % N(i,j,k) sensitivity norm for qoi i, parameter j, day k
    for day = 1:7
        for i = 1:25
            vecplot = zeros(35,500);
            for n = 1:35
                vecplot(n,:) = A(i,:,n,day)';
            end
            for k = 1:35
                sens_norm = sqrt(sum(vecplot.^2,2));
                N(i,:,day) = sens_norm;
            end
        end
    end
end