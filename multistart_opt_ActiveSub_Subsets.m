%% optimizing using the randomized patient data over active variables
clear; clc; close all;
parameterNames = {'Rmval','Raval','Rtval','Rpval','Rs','Rvc','Rp','Rpv' ...
'Csa','Csv','Cpa','Cpv', ...
'Emax_LA','Emin_LA','Vd_LA','T_peak_LA','T_relax_LA','t_offset_LA', ...
'Emax_LV','Emin_LV','Vd_LV','T_peak_LV','T_relax_LV', ...
'Emax_RA','Emin_RA','Vd_RA','T_peak_RA','T_relax_RA','t_offset_RA', ...
'Emax_RV','Emin_RV','Vd_RV','T_peak_RV','T_relax_RV'};
%%
% Specify which patient to analyze
which_patient=2;
% Set the number of initializations for multistart optimization
num_init = 6;
% Decide which parameters you are fixing
param_ids = 1:34;
param_remove = [10 11 12 19]; %Consider all compliance parameters (9-12) and emax in the lV (19)
Vd_remove = [15 21 26 32];
param_ids([param_remove Vd_remove]) = [];

% Determine how many subspaces you want to use
num_act_variables = 2;
% MAKE SURE TO CHANGE THE FILE SAVE NAME
fname = "multioptResults_ActSub_TESTING_pat2_6init";
%%
disp('Not including these parameters in active subspace');
disp(parameterNames(param_remove))
parameterNames(param_remove) = [];

%%
[IC, pars] = ComputeParameters_adult(70,0.0,1);
HBs = 1;

%% Weights: we can specify how much we would like to weigh parts of
% the cost functional
w1 = 1;
w2 = 1;
cost_weights = [w1.*ones(8,1);w2.*ones(6,1)];
%%
load("randomizedPatientData");
load("Samples_for_AS.mat","random_parameter_array")
data = squeeze(patient_group_A(:,:,which_patient));
W = generate_AS_fromData(data,cost_weights,param_ids);



UB = log(2.*max(random_parameter_array,[],2));
LB = log(0.5.*min(random_parameter_array,[],2));

pars_base = makearray(pars);
pars_base([end]) = [];

AS_min = -15;
AS_max = 5;

new_UB = [UB(param_remove);AS_max.*ones(num_act_variables,1)];
new_LB = [LB(param_remove);AS_min.*ones(num_act_variables,1)];

%% active subspaces
% HERE WE ASSUME EVERYTHING IS IN LOG SCALE
pars_base = log(pars_base);
w_active = W(:,1:num_act_variables);
w_inact  = W(:,num_act_variables+1:end);

n_inact    = size(w_inact,2);
n_normparam = length(param_remove);



%% initial guesses
alpha = w_active'*pars_base(param_ids);
beta = w_inact'*pars_base(param_ids);

w_active_full = zeros(34,num_act_variables); 
w_inact_full = zeros(34,n_inact);
w_active_full(param_ids,1:num_act_variables) = w_active;
w_inact_full(param_ids,1:n_inact) = w_inact;

%%
options = optimoptions("fmincon","Display","iter-detailed", ...
    "FiniteDifferenceStepSize",1e-3);
results = zeros(num_act_variables+length(param_remove),num_init);
param_results = zeros(34,num_init); % column i is the inferred parameters for patient i from minimizing cost over the active variables
costs = zeros(1,num_init);

%% Define a function that can take in parameters and active subspace variables


% functions to optimize 1d active subspace and 2d active subspace
% cost_eightCompModel_withSubsets(data,IC,theta_in,param_ids,Wact,Winact,beta,HBs,T,cost_weights,plot_flag)
J = @(x) cost_eightCompModel_withSubsets(data,IC,x, ...
    param_ids,param_remove,w_active,w_inact,beta,...
    HBs,pars.T,cost_weights,0); 
%% Optimization
% Initalize and run multiple optimization
param_init = [pars_base(param_remove);alpha];
for i=1:num_init
    param0 = param_init.*unifrnd(0.75,1.25,length(param_init),1); % Perturb initial value
    disp([i;param0]')
    [xmin,fval,exitflag,output,lambda,grad,hessian] = fmincon(J,param0,[],[],[],[],new_LB,new_UB,[],options);

    results(:,i) = xmin;
    theta_normal = xmin(1:length(param_remove));
    x_AS         = xmin(length(param_remove)+1:end);

    theta_vec = w_active_full*x_AS+w_inact_full*beta;
    theta_vec(param_remove) = theta_normal;

    param_results(:,i) = exp(theta_vec);

    costs(:,i) = fval;
    
    % This part of the script generates a file plot for you to observe the
    % results
    cost_eightCompModel_withSubsets(data,IC,xmin, ...
    param_ids,param_remove,w_active,w_inact,beta,...
    HBs,pars.T,cost_weights,1);
end
    cost_eightCompModel_withSubsets(data,IC,param_init, ...
    param_ids,param_remove,w_active,w_inact,beta,...
    HBs,pars.T,cost_weights,1);

%%
parameterNames = {'Rmval','Raval','Rtval','Rpval','Rs','Rvc','Rp','Rpv' ...
    'Csa','Csv','Cpa','Cpv', ...
    'Emax_LA','Emin_LA','Vd_LA','T_peak_LA','T_relax_LA','t_offset_LA', ...
    'Emax_LV','Emin_LV','Vd_LV','T_peak_LV','T_relax_LV', ...
    'Emax_RA','Emin_RA','Vd_RA','T_peak_RA','T_relax_RA','t_offset_RA', ...
    'Emax_RV','Emin_RV','Vd_RV','T_peak_RV','T_relax_RV'};

param_true = random_patient_parameters(:,which_patient);
rel_par_error = (param_results - param_true)./param_true;
figure(100);
bar(rel_par_error);
xticks(1:34);
xticklabels(parameterNames)
%%
save(fname,"results","param_results",...
    "alpha","beta","w_active","w_inact");