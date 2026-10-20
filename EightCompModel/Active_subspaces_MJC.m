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
num_outputs = 16;
num_param = 34;
num_sample = 3000;
h = 1e-6;

[IC,pars] = ComputeParameters_adult(70,sex);
T = pars.T;
cost_weights = ones(num_outputs,1);
plot_flag = 0;
[V,p,q,tplot] = GetModelOutput(IC,pars,HBs);

% converting the parameters to an array for convenience later
parscell = struct2cell(pars);
pars_array = cell2mat(parscell);

volumescell = struct2cell(V);
pressurescell = struct2cell(p);
flowscell = struct2cell(q);

volumes = cell2mat(volumescell);
pressures = cell2mat(pressurescell);
flows = cell2mat(flowscell);

% 24 by 500 matrix containing the reference subject quantities of interest at 500 time points
data = [volumes; pressures; flows];

%% sampling parameters
% number of samples, step size for forward difference, % width of intervals
% for parameter values
percent_width = 0.50;

L = (1-percent_width)*pars_array(1:num_param,1);
U = (1+percent_width)*pars_array(1:num_param,1);

random_parameter_array = zeros(num_param,num_sample);
cost_evals = zeros(num_sample,1);
model_evals = zeros(num_outputs,num_sample);
grad_evals  = zeros(num_outputs,num_param,num_sample);
for i = 1:num_sample
    random_parameter_array(:,i) = unifrnd(L,U,[num_param,1]);
end


%% computing C
% empty monte carlo approximation
C = zeros(num_param,num_param);

tic
parfor i = 1:num_sample
    disp(i)
    grad = zeros(num_param,1);
    grad_qoi  = zeros(num_outputs,num_param);
    theta = pars;
    rand_pars = pars;
    % put random parameters into a struct to pass to the model
    % (not sampling a random VSD width)
    rand_pars_vector = log(random_parameter_array(:,i));

    % rand_Rpars = [rand_pars_vector(1:8)];
    % rand_Cpars = [rand_pars_vector(9:12)];
    % rand_LApars = [rand_pars_vector(13:18)];
    % rand_LVpars = [rand_pars_vector(19:23)];
    % rand_RApars = [rand_pars_vector(24:29)];
    % rand_RVpars = [rand_pars_vector(30:34)];
    % 
    % rand_pars.Rpars = rand_Rpars;
    % rand_pars.Cpars = rand_Cpars;
    % rand_pars.LApars = rand_LApars;
    % rand_pars.LVpars = rand_LVpars;
    % rand_pars.RApars = rand_RApars;
    % rand_pars.RVpars = rand_RVpars;
    % rand_pars.T = pars.T;

    % cost function value at parameter sample i
    [ref_cost,ref_qoi] = cost_eightCompModel_MJC(data,IC,rand_pars_vector,HBs,T,cost_weights,plot_flag);
    cost_evals(i)    = ref_cost;
    model_evals(:,i) = ref_qoi;

    % finite difference in each direction for 35 parameters
    for j = 1:num_param
        perturbed_pars_vector = rand_pars_vector;

        % perturbed_pars_array = pars_array*1; % unperturb parameter j-1
        perturbed_pars_vector(j) = perturbed_pars_vector(j) + h; % perturb parameter j

        % put the parameters back into structure form
        % rand_Rpars = [perturbed_pars_vector(1:8)];
        % rand_Cpars = [perturbed_pars_vector(9:12)];
        % rand_LApars = [perturbed_pars_vector(13:18)];
        % rand_LVpars = [perturbed_pars_vector(19:23)];
        % rand_RApars = [perturbed_pars_vector(24:29)];
        % rand_RVpars = [perturbed_pars_vector(30:34)];
        % theta.Rpars = rand_Rpars;
        % theta.Cpars = rand_Cpars;
        % theta.LApars = rand_LApars;
        % theta.LVpars = rand_LVpars;
        % theta.RApars = rand_RApars;
        % theta.RVpars = rand_RVpars;
        % 
        % theta.T = pars.T;

        % update jth gradient component using forward finite difference
        [new_cost,new_qoi] = cost_eightCompModel_MJC(data,IC,perturbed_pars_vector,HBs,T,cost_weights,plot_flag);
        grad(j) = (new_cost - ref_cost)/h;
        grad_qoi(:,j)  = (new_qoi-ref_qoi)/h;
    end

    % add outer product to monte carlo approximation
    grad_evals(:,:,i) = grad_qoi;
    C = C + grad*grad';

end
toc

% divide by the number of samples
C = 1/num_sample*C;


%% active subspace
% eigendecomposition of C, W eigenvectors, w eigenevalues
[W, w] = eig(C);

save("AS_8Comp_logparam_10_14_26","W","w","C","random_parameter_array","grad_evals",'cost_evals','model_evals');

%% Plotting
