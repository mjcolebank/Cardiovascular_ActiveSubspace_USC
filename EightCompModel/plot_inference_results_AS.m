%% optimizing using the randomized patient data over active variables
clear; clc; close all;
parameterNames = {'Rmval','Raval','Rtval','Rpval','Rs','Rvc','Rp','Rpv' ...
    'Csa','Csv','Cpa','Cpv', ...
    'Emax_LA','Emin_LA','Vd_LA','T_peak_LA','T_relax_LA','t_offset_LA', ...
    'Emax_LV','Emin_LV','Vd_LV','T_peak_LV','T_relax_LV', ...
    'Emax_RA','Emin_RA','Vd_RA','T_peak_RA','T_relax_RA','t_offset_RA', ...
    'Emax_RV','Emin_RV','Vd_RV','T_peak_RV','T_relax_RV'};
%%
% parpool(6);
% Specify which patient to analyze
for which_patient=5%1:5
    for num_act_variables=1%:-1:1
        % Set the number of initializations for multistart optimization
        num_init = 6;
        % Decide which parameters you are fixing
        param_ids = 1:34;
        param_remove = [10 11];%[10 11 12 19]; %Consider all compliance parameters (9-12) and emax in the lV (19)
        Vd_remove = [15 21 26 32];
        param_ids([param_remove Vd_remove]) = [];

        % Determine how many subspaces you want to use
        % num_act_variables = 6;
        % MAKE SURE TO CHANGE THE FILE SAVE NAME
        fname = strcat("multiopt_AS_pat",num2str(which_patient),"_NumAS_",num2str(num_act_variables),"_NumTruePar_",num2str(length(param_remove)));
        %%
        disp('Not including these parameters in active subspace');
        disp(parameterNames(param_remove))
        parameterNames(param_remove) = [];

        %%
        [IC, pars] = ComputeParameters_adult(70,0.0,1);
        T = pars.T;
        HBs = 1;

        %% Weights: we can specify how much we would like to weigh parts of
        % the cost functional
        w1 = 1;
        w2 = 1;
        cost_weights = [w1.*ones(8,1);w2.*ones(6,1)];
        %%
        load("randomizedPatientData");
        load("Samples_for_AS_NewTiming.mat","random_parameter_array")
        data = squeeze(patient_group_A(:,:,which_patient));
        W = generate_AS_fromData(data,cost_weights,param_ids);
        %%
        qdata = zeros(14,1);
        % maximum and minimum pressure in the systemic arteries
        qdata(1) = max(data(13,:));qdata(2) = min(data(13,:));

        % maximum and minimum volume in the left ventricle
        qdata(3) = max(data(2,:));qdata(4) = min(data(2,:));

        % maximum and minimum volume in the left atrium
        qdata(5) = max(data(1,:));qdata(6) = min(data(1,:));

        % maximum and minimum volume in the right ventricle
        qdata(7) = max(data(6,:));qdata(8) = min(data(6,:));


        % maximum flow rates in the four valves mitral, aortic, tricuspid,
        % pulmonary
        qdata(9) = max(data(17,:));
        qdata(10) = max(data(18,:));
        qdata(11) = max(data(19,:));
        qdata(12) = max(data(20,:));

        % mean right atrial pressure
        qdata(13) = mean(data(11,:));
        % maximum pulmonary artery pressure
        qdata(14) = max(data(15,:));
        %%



        UB = log(2.*max(random_parameter_array,[],2));
        LB = log(0.5.*min(random_parameter_array,[],2));

        pars_base = makearray(pars);
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
        % alpha = w_active'*pars_base(param_ids);
        % beta = w_inact'*pars_base(param_ids);
        load(fname)
    
        w_active_full = zeros(34,num_act_variables);
        w_inact_full = zeros(34,n_inact);
        w_active_full(param_ids,1:num_act_variables) = w_active;
        w_inact_full(param_ids,1:n_inact) = w_inact;
        

        
        %% Optimization
        % Initalize and run multiple optimization
        param_init = [pars_base(param_remove);alpha];
        for i=1:num_init
            param_curr = makestruct(param_results(:,i),10e8,T);
            [V,p,q,tplot] = GetModelOutput(IC,param_curr,HBs);
            qoi = zeros(14,1);




            figure(100+i);
            subplot(2,2,1); hold on;
            plot(data(1,:),data(9,:),'k');
            plot(V.LA,p.LA,'r');
            xline(qdata(5:6),'--k');

            subplot(2,2,2); hold on;
            plot(data(2,:),data(10,:),'k');
            plot(V.LV,p.LV,'r');
            xline(qdata(3:4),'--k');

            subplot(2,2,3); hold on;
            plot(p.SA,'r');
            plot(p.RA,'c')
            plot(p.PA,'b');
            yline(qdata(1:2),'--r');
            yline(qdata(13),'--c');
            yline(qdata(14),'--b');

            subplot(2,2,4); hold on;
            plot(data(6,:),data(12,:),'k');
            plot(V.RV,p.RV,'b');
            xline(qdata(7:8),'--k');

            figure(200+i);
            subplot(2,2,1); hold on;
            plot(data(17,:),'k');
            plot(q.MVal,'r');
            yline(qdata(9),'--k')

            subplot(2,2,2); hold on;
            plot(data(18,:),'k');
            plot(q.AVal,'r');
            yline(qdata(10),'--k')

            subplot(2,2,3); hold on;
            plot(data(19,:),'k');
            plot(q.TVal,'r');
            yline(qdata(11),'--k')

            subplot(2,2,4); hold on;
            plot(data(20,:),'k');
            plot(q.PVal,'r');
            yline(qdata(12),'--k')



        end
        

        %%
        parameterNames = {'Rmval','Raval','Rtval','Rpval','Rs','Rvc','Rp','Rpv' ...
            'Csa','Csv','Cpa','Cpv', ...
            'Emax_LA','Emin_LA','Vd_LA','T_peak_LA','T_relax_LA','t_offset_LA', ...
            'Emax_LV','Emin_LV','Vd_LV','T_peak_LV','T_relax_LV', ...
            'Emax_RA','Emin_RA','Vd_RA','T_peak_RA','T_relax_RA','t_offset_RA', ...
            'Emax_RV','Emin_RV','Vd_RV','T_peak_RV','T_relax_RV'};

        param_true = random_patient_parameters(:,which_patient);
        param_true([10 end]) = [];
        rel_par_error = (param_results - param_true)./param_true;
        figure(10000);
        bar(rel_par_error);
        xticks(1:34);
        xticklabels(parameterNames)

    end
end