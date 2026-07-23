% Optimization within active subspace
% Author: Ethan Mahan
% 6/14 2026

%% reference data

VSD_rad = 0.0; % cm
HBs = 1;
age = 7; % days
sex = 1;
percentile = 50;

[IC,pars] = ComputeParameters_adult(70,VSD_rad,sex);

[V,p,q,plot] = GetModelOutput(IC,pars,HBs);

% put the parameters in a matrix for convenience later
parscell = struct2cell(pars);
pars_array = cell2mat(parscell); 

volumescell = struct2cell(V);
pressurescell = struct2cell(p);
flowscell = struct2cell(q);

volumes = cell2mat(volumescell);
pressures = cell2mat(pressurescell);
flows = cell2mat(flowscell);

% 25 by 500 matrix containing the reference quantities of interest at 500 time points
data = [volumes; pressures; flows];

%% optimization setup
% active subspace
load("active_subspace_eightComp_1000samp.mat")
w_active = [W(1:8,34); 0; W(9:34,34)];

% function to optimize
J = @(x) cost2_eightCompModel(data,IC,w_active*x,HBs,pars.T);

% initial guess
y = w_active'*pars_array(1:35);

options = optimoptions("fmincon","Display","iter-detailed");
[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(J,y,[],[],[],[],-1000,1000,[],options);

par_vector = w_active*x;

%% plotting output result
par_struct = makestruct(par_vector,1e10,pars.T);

[V,p,q,tplot] = GetModelOutput(IC,par_struct,HBs);