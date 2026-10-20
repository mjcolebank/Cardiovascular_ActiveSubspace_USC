function parameter_structure = makestruct(A,T)
% takes in 34 d parameter vector and makes it into a structure
% fixed value for VSD_resistance and time cycle T

    parameter_structure.Rpars = [A(1:8)];
    parameter_structure.Cpars = [A(9:12)];
    parameter_structure.LApars = [A(13:18)];
    parameter_structure.LVpars = [A(19:23)];
    parameter_structure.RApars = [A(24:29)];
    parameter_structure.RVpars = [A(30:34)];
    parameter_structure.T = T;

end