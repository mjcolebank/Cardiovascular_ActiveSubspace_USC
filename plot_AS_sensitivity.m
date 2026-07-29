%% optimizing using the randomized patient data over active variables
clear; clc; %close all;
which_patient=10;
num_init = 10;
parameterNames = {'Rmval','Raval','Rtval','Rpval','Rs','Rvc','Rp','Rpv' ...
'Csa','Csv','Cpa','Cpv', ...
'Emax_LA','Emin_LA','Vd_LA','T_peak_LA','T_relax_LA','t_offset_LA', ...
'Emax_LV','Emin_LV','Vd_LV','T_peak_LV','T_relax_LV', ...
'Emax_RA','Emin_RA','Vd_RA','T_peak_RA','T_relax_RA','t_offset_RA', ...
'Emax_RV','Emin_RV','Vd_RV','T_peak_RV','T_relax_RV'};
%%
[IC, pars] = ComputeParameters_adult(70,0.0,1);
HBs = 1;

% load("active_subspace_eightComp_1000samp.mat");
% load("AS_MJC.mat")
load("randomizedPatientData");
load("AS_MJC.mat","random_parameter_array")
data = squeeze(patient_group_A(:,:,which_patient));
[W,lambda] = generate_AS_fromData(data);

% Choose active subspace dimension k
k = length(lambda);   % example
% Compute activity scores
activity_scores = sum(W(:,1:k).^2 .* lambda(1:k)', 2);

% figure; bar(real(log(activity_scores)));
figure; bar(activity_scores);

xticks(1:length(lambda));
xticklabels(parameterNames);
grid on
xtickangle(45);