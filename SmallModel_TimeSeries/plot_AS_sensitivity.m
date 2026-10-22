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
w2 = 1./sqrt(50);%1;
cost_weights = [w1.*ones(2,1);w2.*ones(50,1)];
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
    disp([W(where(1:5),i)]')

end