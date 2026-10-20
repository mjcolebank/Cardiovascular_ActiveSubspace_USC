%% optimizing using the randomized patient data over active variables
clear; clc; close all;
% Specify which patient to analyze
which_patient=1;
% Set the number of initializations for multistart optimization
num_init = 1;
% Determine how many subspaces you want to use
num_act_variables = 4;
% MAKE SURE TO CHANGE THE FILE SAVE NAME
fname = "multioptResults_ActSub_TESTING";
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
W = generate_AS_fromData(data,cost_weights);



UB = log(max(random_parameter_array,[],2));
LB = log(min(random_parameter_array,[],2));

pars_base = makearray(pars);
pars_base([end]) = [];

AS_min = -15;
AS_max = 5;


%% active subspaces
% HERE WE ASSUME EVERYTHING IS IN LOG SCALE
pars_base = log(pars_base);
w_active = W(:,1:num_act_variables);
w_inact  = W(:,num_act_variables+1:end);


%% initial guesses
alpha = w_active'*pars_base;
beta = w_inact'*pars_base;

%%
options = optimoptions("fmincon","Display","iter-detailed", ...
    "FiniteDifferenceStepSize",1e-2);
results = zeros(num_act_variables,num_init);
param_results = zeros(34,num_init); % column i is the inferred parameters for patient i from minimizing cost over the active variables
costs = zeros(1,num_init);

% functions to optimize 1d active subspace and 2d active subspace
J = @(x) cost_eightCompModel_MJC(data,IC,w_active*x+w_inact*beta, ...
    HBs,pars.T,cost_weights,0); 
%% Optimization
% Initalize and run multiple optimization
AS_LB = AS_min.*ones(num_act_variables,1);
AS_UB = AS_max.*ones(num_act_variables,1);

for i=1:num_init
    alpha0 = alpha.*unifrnd(0.5,1.5,num_act_variables,1); % Perturb initial value
    if any(alpha0>AS_max)
        id = find(alpha0>AS_max);
        alpha0(id) = AS_max*0.999;
        disp('hit max')
    end
    if any(alpha0<AS_min)
        id = find(alpha0>AS_max);
        alpha0(id) = AS_min*0.999;
        disp('hit min')
    end
    disp([i;alpha0]')
    [xmin,fval,exitflag,output,lambda,grad,hessian] = fmincon(J,alpha0,[],[],[],[],AS_LB,AS_UB,[],options);

    results(:,i) = xmin;

    param_results(:,i) = exp(w_active*xmin+w_inact*beta);

    costs(:,i) = fval;
    
    % This part of the script generates a file plot for you to observe the
    % results
    cost_eightCompModel_MJC(data,IC,w_active*xmin+w_inact*beta, ...
        HBs,pars.T,cost_weights,1); % x is a scalar
end
    cost_eightCompModel_MJC(data,IC,w_active*alpha+w_inact*beta, ...
        HBs,pars.T,cost_weights,1); % x is a scalar

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