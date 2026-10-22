%% Try to derive the active subspace confidence intervals
clear; clc; close all;
%%
which_patient=5;
how_many_eigs = 5;
parameterNames = {'P_{LA}','P_{SysCap}',...
    'R_{MVal}','R_{AVal}','R_{art}', ...
    'E_{es}','A','B','T_{peak}','T_{relax}',...
    'V_{LV,d}','C_{Ao}'};
n_pars = length(parameterNames);
param_remove = [];
n_outputs = 6;
%%
[IC, pars] = ComputeParameters_adult();
HBs = 1.0;
T = pars(end);
%%
% parpool(6);
% Specify which patient to analyze
for which_patient=1
    for num_act_variables=3%1:3

        % MAKE SURE TO CHANGE THE FILE SAVE NAME
        fname = strcat("multiopt_AS_pat",num2str(which_patient),"_NumAS_",num2str(num_act_variables),"_NumTruePar_",num2str(length(param_remove)));
        % fname = strcat("multiopt_noAS_pat",num2str(which_patient));
        load(fname);
        which_opt = 1;
        param_opt = results(:,which_opt);
        alpha_opt = alpha;
        beta_opt = beta;
        Wa = w_active;
        Wb = w_inact;
        %% Weights: we can specify how much we would like to weigh parts of
        % the cost functional
        w1 = 1;
        w2 = 1;
        w3 = 1;
        cost_weights = [w1.*ones(2,1);w2.*ones(2,1);w3.*ones(2,1)];
        %%
        disp('No longer including these parameters');
        disp(parameterNames(param_remove))
        parameterNames(param_remove) = [];

        %%

        load("randomizedPatientData");
        data = squeeze(patient_group_A(:,:,which_patient));
        [W,lambda] = generate_AS_fromData(data,cost_weights,param_ids);
        
        h = 1e-6;
        n_param = length(param_opt);
        grad_qoi  = zeros(n_outputs,n_param);
        [ref_cost,ref_qoi] = cost_SmallModel_NOSCALE(data,IC,param_opt,param_ids,param_remove,Wa,Wb,beta,HBs,T,1);
        for k = 1:n_param
            perturbed_pars_array = param_opt;
            perturbed_pars_array(k) = perturbed_pars_array(k) + h; % perturb parameter j
            % update jth gradient component using forward finite difference
            % theta_eval = Wa*exp(perturbed_pars_array) + Wb*beta;
            [new_cost,new_qoi] = cost_SmallModel_NOSCALE(data,IC,perturbed_pars_array,param_ids,param_remove,Wa,Wb,beta,HBs,T,0);            
            grad_qoi(:,k)  = (new_qoi-ref_qoi)/h;
            
            % grad(j) = (log(new_cost) - log(ref_cost))/h;
            % grad_qoi(:,j)  = (log(new_qoi)-log(ref_qoi))/h;
        end

        % Now for confidence interval calculations
        s2_est = 1./(n_outputs-n_param).*ref_cost;

    end
end