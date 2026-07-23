%% optimizing using the randomized patient data over active variables
clear; clc; close all;

%%
[IC, pars] = ComputeParameters_adult(70,0.0,1);
HBs = 1;

load("active_subspace_eightComp_1000samp.mat");
load("randomizedPatientData");

UB = max(random_parameter_array,[],2);
LB = min(random_parameter_array,[],2);

pars_base = makearray(pars);
pars_base(end) = [];

%% active subspaces
W = [W(1:8,:); zeros(1,34); W(9:34,:)];
W = [W(:,1:8), zeros(35,1), W(:,9:34)];

w_active_1D = W(:,end);
w_inact_1D  = W(:,1:end-1);

w_active_2D = W(:,end-1:end);
w_inact_2D  = W(:,1:end-2);

w_active_3D = W(:,end-2:end);
w_inact_3D = W(:,1:end-3);

%% initial guesses
alpha_1D = w_active_1D'*pars_base;
beta_1D = w_inact_1D'*pars_base;
list = linspace(alpha_1D*0.95,alpha_1D*1.05,20);

alpha_2D = w_active_2D'*pars_base;
beta_2D = w_inact_2D'*pars_base;
x2_list = linspace(alpha_2D(1)*0.95,alpha_2D(1)*1.05,20);
y2_list = linspace(alpha_2D(2)*0.95,alpha_2D(2)*1.05,20);

alpha_3D = w_active_3D'*pars_base;
beta_3D = w_inact_3D'*pars_base;
x3_list = linspace(alpha_2D(1)*0.95,alpha_2D(1)*1.05,20);
y3_list = linspace(alpha_2D(2)*0.95,alpha_2D(2)*1.05,20);
z3_list = linspace(alpha_3D(3)*0.95,alpha_3D(3)*1.05,20);

vsd_offset = zeros(35,1);
vsd_offset(9) =  1e10;

%%
% w_active_10D = W(:,25:35);
options = optimoptions("fmincon","Display","iter-detailed");
results1D = zeros(number_of_patients,1);
results2D = zeros(number_of_patients,2); % each row is the 2d vector result from optimizer for the ith patient
param_results1D = zeros(35,number_of_patients); % column i is the inferred parameters for patient i from minimizing cost over the active variables
param_results2D = zeros(35,number_of_patients);
costs1D = zeros(number_of_patients,1);
costs2D = zeros(number_of_patients,1);

% patient data, cost function reference to each patient i
patient_data = squeeze(patient_group_A(:,:,1));

% functions to optimize 1d active subspace and 2d active subspace
J_1 = @(x) cost2_eightCompModel(patient_data,IC,w_active_1D*x+w_inact_1D*beta_1D+vsd_offset,HBs,pars.T); % x is a scalar
J_2 = @(x2) cost2_eightCompModel(patient_data,IC,w_active_2D*x2+w_inact_2D*beta_2D+vsd_offset,HBs,pars.T); % x2 is a 2d vector
J_3 = @(x) cost2_eightCompModel(patient_data,IC,w_active_3D*x+w_inact_3D*beta_3D+vsd_offset,HBs,pars.T);

% polling
answers_J2 = zeros(20,20);
answers_J3 = zeros(20,20,20);

for i = 1:20
    for j=1:20
        answers_J2(i,j) = J_2([x2_list(i);y2_list(j)]);
        for k=1:20
            answers_J3(i,j,k) = J_3([x3_list(i);y3_list(j);z3_list(k)]);
        end
    end
end

%%
% results
[xmin,fval,exitflag,output,lambda,grad,hessian] = fmincon(J_1,alpha_1D,[],[],[],[],20,40,[],options);
[x2min,fval2,exitflag2,output2,lambda2,grad2,hessian2] = fmincon(J_2,alpha_2D,[],[],[],[],[-1000,-1000],[1000,1000],[],options);

results1D(i) = xmin;
% results2D(i,:) = x2min;

param_results1D(:,i) = w_active_1D*xmin;
% param_results2D(:,i) = w_active_2D*x2min;

costs1D(i) = fval;
% costs2D(i) = fval2;

save("optimizationResults","results1D","results2D","param_results1D","param_results2D");