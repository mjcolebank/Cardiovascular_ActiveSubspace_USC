function parameter_structure = makestruct(A,VSD_resistance,T)
% takes in 35 d parameter vector and makes it into a structure
% fixed value for VSD_resistance and time cycle T

    parameter_structure.Rpars = [A(1:8); VSD_resistance];
    parameter_structure.Cpars = [A(10:13)];
    parameter_structure.LApars = [A(14:19)];
    parameter_structure.LVpars = [A(20:24)];
    parameter_structure.RApars = [A(25:30)];
    parameter_structure.RVpars = [A(31:35)];
    parameter_structure.T = T;

end