%% optimizing using the randomized patient data over active variables
clear; clc; close all;
which_patient=2;
how_many_eigs = 4;
parameterNames = {'Rmval','Raval','Rtval','Rpval','Rs','Rvc','Rp','Rpv' ...
    'Csa','Csv','Cpa','Cpv', ...
    'Emax_LA','Emin_LA','Vd_LA','T_peak_LA','T_relax_LA','t_offset_LA', ...
    'Emax_LV','Emin_LV','Vd_LV','T_peak_LV','T_relax_LV', ...
    'Emax_RA','Emin_RA','Vd_RA','T_peak_RA','T_relax_RA','t_offset_RA', ...
    'Emax_RV','Emin_RV','Vd_RV','T_peak_RV','T_relax_RV'};
%%
[IC, pars] = ComputeParameters_adult(70,1);
HBs = 1;


%% Weights: we can specify how much we would like to weigh parts of
% the cost functional
w1 = 1;
w2 = 1;
cost_weights = [w1.*ones(10,1);w2.*ones(6,1)];
%%
param_ids = 1:34;
param_remove = [];%[9:12];%[10 11 12 19]; %Consider all compliance parameters (9-12) and emax in the lV (19)
% param_remove = param_ids; param_remove(9:12) = [];
% Also add dead volumes
param_remove(end+1:end+4) = [15 21 26 32];
param_ids(param_remove) = [];
disp('No longer including these parameters');
disp(parameterNames(param_remove))
parameterNames(param_remove) = [];

%%

load("randomizedPatientData");
data = squeeze(patient_group_A(:,:,which_patient));
[W,lambda] = generate_AS_fromData(data,cost_weights,param_ids);

%%
qdata = zeros(16,1);
% maximum and minimum pressure in the systemic arteries
qdata(1) = max(data(13,:));
qdata(2) = min(data(13,:));
% maximum and minimum volume in the left ventricle
qdata(3) = max(data(2,:));
qdata(4) = min(data(2,:));
% maximum and minimum volume in the left atrium
qdata(5) = max(data(1,:));
qdata(6) = min(data(1,:)); 
% maximum and minimum volume in the right ventricle
qdata(7) = max(data(6,:));
qdata(8) = min(data(6,:));
% maximum and minimum volume in the right atrium
qdata(9) = max(data(5,:));
qdata(10) = min(data(5,:));

% maximum flow rates in the four valves mitral, aortic, tricuspid,
% pulmonary
qdata(11) = max(data(17,:));
qdata(12) = max(data(18,:));
qdata(13) = max(data(19,:));
qdata(14) = max(data(20,:));

% mean right atrial pressure
qdata(15) = mean(data(11,:));

% maximum pulmonary artery pressure
qdata(16) = max(data(15,:));

%%
load AS_8Comp_logparam_10_14_26.mat model_evals  random_parameter_array cost_evals
random_parameter_array(param_remove,:) = [];
J = sum(((qdata - model_evals)./qdata).^2);

%% active subspaces
% HERE WE ASSUME EVERYTHING IS IN LOG SCALE
pars_base = log(random_parameter_array);
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


