%% optimizing using the randomized patient data over active variables
clear; clc; %close all;
which_patient=5;
how_many_eigs = 5;
parameterNames = {'P_{LA}','P_{SysCap}',...
    'R_{MVal}','R_{AVal}','R_{art}', ...
    'E_{es}','A','B','T_{peak}','T_{relax}',...
    'V_{LV,d}','C_{Ao}'};
%%
[IC, pars] = ComputeParameters_adult();
HBs = 1;


%% Weights: we can specify how much we would like to weigh parts of
% the cost functional
w1 = 1;
w2 = 1;
w3 = 1;
cost_weights = [w1.*ones(2,1);w2.*ones(2,1);w3.*ones(2,1)];
%%
param_ids = 1:12;
param_remove = [];%Consider all compliance parameters (9-12) and emax in the lV (19)
% Also add dead volumes
param_ids(param_remove) = [];

disp('No longer including these parameters');
disp(parameterNames(param_remove))
parameterNames(param_remove) = [];

%%

load("randomizedPatientData");
data = squeeze(patient_group_A(:,:,which_patient));
[W,lambda] = generate_AS_fromData(data,cost_weights,param_ids);

%%
qdata = zeros(6,1);
% maximum and minimum pressure in the systemic arteries
qdata(1) = max(data(4,:));
qdata(2) = min(data(4,:));
% maximum and minimum volume in the left ventricle
qdata(3) = max(data(1,:));
qdata(4) = min(data(1,:));
% maximum flow in the MV and AV valves
qdata(5) = max(data(5,:));
qdata(6) = max(data(6,:)); 


%%
load AS_SmallModel_log_ODE23s_10_13_26.mat model_evals  random_parameter_array cost_evals
random_parameter_array(param_remove,:) = [];
J = sum(((qdata - model_evals)./qdata).^2);

%% active subspaces
% HERE WE ASSUME EVERYTHING IS IN LOG SCALE
pars_base = random_parameter_array;%log(random_parameter_array);
w_active = W(:,1:5);
w_inact  = W(:,5+1:end);
alpha = w_active'*pars_base;
beta = w_inact'*pars_base;
for i=1:5
    % figure(i);clf;
    % for j=1:14
    %     subplot(3,5,j);
    %     plot(alpha(i,:),model_evals(j,:),'o');
    % end
    figure(100+i);clf;
    plot(alpha(i,:),log(J),'o');
end


