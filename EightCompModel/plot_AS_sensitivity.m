%% optimizing using the randomized patient data over active variables
clear; clc; %close all;
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
param_remove = [9 10];%[9:12];%[9:12];%[10 11 12 19]; %Consider all compliance parameters (9-12) and emax in the lV (19)
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

% Choose active subspace dimension k
k = length(lambda);   % example
% Compute activity scores
activity_scores_unnorm = sum(W(:,1:k).^2 .* lambda(1:k)', 2);
activity_scores = activity_scores_unnorm./sum(lambda);
[who,where] = sort(activity_scores,'descend');


figure; bar(activity_scores(where));
xticks(1:length(lambda));
xticklabels(parameterNames(where));
grid on
xtickangle(45);
set(gca,'FontSize',20);

%% Plot eigenvectors corresponding to eigenvalues
% Want to verify that the dominant eigenvectors are not unit vectors
figure;
bar(W(:,1:how_many_eigs));
xticks(1:length(lambda));
xticklabels(parameterNames);
grid on
xtickangle(45);
set(gca,'FontSize',20);
legend('W1','W2','W3')
for i=1:how_many_eigs
    [who,where] = sort(abs(W(:,i)),'descend');
    disp('-------------------------------------')
    disp([parameterNames(where(1:5))])
    disp([who(1:5)]')

end