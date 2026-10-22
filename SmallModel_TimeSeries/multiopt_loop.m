%% optimizing using the randomized patient data over active variables
clear; clc; close all;
parameterNames = {'P_{LA}','P_{SysCap}',...
    'R_{MVal}','R_{AVal}','R_{art}', ...
    'E_{es}','A','B','T_{peak}','T_{relax}',...
    'V_{LV,d}','C_{Ao}'};
n_pars = length(parameterNames);
%%
 % parpool(6); 
% Specify which patient to analyze
for which_patient=1
    for num_act_variables=1:6
        % Set the number of initializations for multistart optimization
        num_init = 4;
        % Decide which parameters you are fixing
        param_ids = 1:n_pars;
        param_remove = [];%[10 11 12 19]; %Consider all compliance parameters (9-12) and emax in the lV (19)
        param_ids(param_remove) = [];

        % Determine how many subspaces you want to use
        % num_act_variables = 6;
        % MAKE SURE TO CHANGE THE FILE SAVE NAME
        fname = strcat("multiopt_AoTimeSeries_AS_pat",num2str(which_patient),"_NumAS_",num2str(num_act_variables),"_NumTruePar_",num2str(length(param_remove)));
        % fname = strcat("multiopt_noAS_pat",num2str(which_patient));

        %%
        disp('Not including these parameters in active subspace');
        disp(parameterNames(param_remove))
        parameterNames(param_remove) = [];

        %%
        [IC, pars] = ComputeParameters_adult();
        HBs = 1;
        T = pars(end);
        %% Weights: we can specify how much we would like to weigh parts of
        % the cost functional
        w1 = 1;%1./93;
        w2 = 1;%1./100;
        w3 = 1;%1./400;
        cost_weights = [w1.*ones(2,1);w2.*ones(50,1)];
        %%
        load("randomizedPatientData");
        load("AS_AoTimeSeries_SmallModel_log_ODE23s_10_13_26.mat","random_parameter_array")
        data = squeeze(patient_group_A(:,:,which_patient));
        W = generate_AS_fromData(data,cost_weights,param_ids);



        UB = log(2.*max(random_parameter_array,[],2));
        LB = log(0.5.*min(random_parameter_array,[],2));

        pars_base = pars';
        pars_base([end]) = [];

        AS_min = -15;
        AS_max = 15;

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

        w_active_full = zeros(n_pars,num_act_variables);
        w_inact_full = zeros(n_pars,n_inact);
        w_active_full(param_ids,1:num_act_variables) = w_active;
        w_inact_full(param_ids,1:n_inact) = w_inact;

        %%
        opt_opt = optimoptions("fmincon","Display","iter-detailed", ...
            "FiniteDifferenceStepSize",1e-3,"UseParallel",true);
        results = zeros(num_act_variables+length(param_remove),num_init);
        param_results = zeros(n_pars,num_init); % column i is the inferred parameters for patient i from minimizing cost over the active variables
        costs = zeros(1,num_init);

        %% Define a function that can take in parameters and active subspace variables


        % functions to optimize 1d active subspace and 2d active subspace
        % cost_eightCompModel_withSubsets(data,IC,theta_in,param_ids,Wact,Winact,beta,HBs,T,cost_weights,plot_flag)
        J = @(x) cost_SmallCompModel_withSubsets(data,IC,x, ...
            param_ids,param_remove,w_active,w_inact,beta,...
            HBs,T,cost_weights,0);
        %% Optimization
        % Initalize and run multiple optimization
        param_init = [pars_base(param_remove);alpha];
        for i=1:num_init
            param0 = param_init.*unifrnd(0.75,1.25,length(param_init),1); % Perturb initial value
            disp([i;param0]')
            [xmin,fval,exitflag,output,lambda,grad,hessian] = fmincon(J,param0,[],[],[],[],new_LB,new_UB,[],opt_opt);

            results(:,i) = xmin;
            theta_normal = xmin(1:length(param_remove));
            x_AS         = xmin(length(param_remove)+1:end);

            theta_vec = w_active_full*x_AS+w_inact_full*beta;
            theta_vec(param_remove) = theta_normal;

            param_results(:,i) = exp(theta_vec);

            costs(:,i) = fval;

            % This part of the script generates a file plot for you to observe the
            % results
            cost_SmallCompModel_withSubsets(data,IC,xmin, ...
                param_ids,param_remove,w_active,w_inact,beta,...
                HBs,T,cost_weights,1);
        end
        % cost_eightCompModel_withSubsets(data,IC,param_init, ...
        %     param_ids,param_remove,w_active,w_inact,beta,...
        %     HBs,pars.T,cost_weights,1);

       
        save(fname,"results","param_results",...
            "alpha","beta","w_active","w_inact","param_remove","param_ids");
    end
end