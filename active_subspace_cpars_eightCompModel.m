% Compute active subspace for 8 compartment model
% Author: Ethan Mahan
% 5/18/2026

%% reference data
% data based on a 50th percentile adult male, no VSD

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
  
%% sampling parameters
% number of samples, step size for forward difference, % width of intervals
% for parameter values
m = 1000;
h = 1e-8;
percent_width = 0.2;

% upper and lower bounds for sampling compliance parameters
L = (1-percent_width)*pars_array(10:13,1); 
U = (1+percent_width)*pars_array(10:13,1);
random_parameter_array = zeros(4,m);

% m uniform random samples for compliance parameters
for i = 1:m
    random_parameter_array(:,i) = unifrnd(L,U,[4,1]);
end

%% computing C
C = zeros(4,4);

tic
parfor i = 1:m
    grad = zeros(4,1);

    % theta is going to be the perturbed parameter structure
    theta = pars;
    rand_pars = pars;

    % random sample for the compliance parameters
    rand_pars.Cpars = random_parameter_array(:,i);

    % cost function value at parameter sample i
    ref_cost = cost_eightCompModel(data,IC,rand_pars,HBs);
    
    % finite difference in each direction for 4 compliance parameters
    for j = 1:4
        % unperturbed parameters
        perturbed_Cpars = rand_pars.Cpars; 

        % perturbing parameter j on a log scale
        perturbed_Cpars(j) = exp(log(perturbed_Cpars(j)) + h);

        % perturbed parameters in structure to pass to model
        theta.Cpars = perturbed_Cpars;

        % finite difference for compliance parameter j
        grad(j) = (cost_eightCompModel(data,IC,theta,HBs) - ref_cost)/h;
    end

    % add outer product to monte carlo approximation
    C = C + grad*grad';
    
end
toc

% divide by the number of samples
C = 1/m*C;

%% active subspace
% eigendecomposition of C, W eigenvectors, w eigenevalues
[W, w] = eig(C);

% sum(W.^2*diag(w),2) vector with activity scores