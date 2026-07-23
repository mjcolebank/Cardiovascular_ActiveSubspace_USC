%% optimizing using the randomized patient data over active variables

load("active_subspace_eightComp_1000samp.mat");
load("randomizedPatientData");

w_active_1D = [W(1:8,34); 0; W(9:34,34)];
w_active_2D = [W(1:8,33), W(1:8,34); 
               0, 0;
               W(9:34,33), W(9:34,34)];

W = [W(1:8,:); zeros(1,34); W(9:34,:)];
W = [W(:,1:8), zeros(35,1), W(:,9:34)];


%%
w_active_10D = W(:,25:35);

% functions to optimize 1d active subspace and 2d active subspace
% J_1 = @(x) cost2_eightCompModel(data,IC,w_active_1D*x,HBs,pars.T); % x is a scalar
% J_2 = @(x2) cost2_eightCompModel(data,IC,w_active_2D*x2,HBs,pars.T); % x2 is a 2d vector

% initial guesses
% y_1 = w_active_1D'*random_parameter_array(1:35,i);
% y_2 = w_active_2D;

% y = 0.1;
y2 = [100; 100;];

% guess_array_1D = [];
% guess_array_2D = [];



options = optimoptions("fmincon","Display","iter-detailed");
results1D = zeros(number_of_patients,1);
results2D = zeros(number_of_patients,2); % each row is the 2d vector result from optimizer for the ith patient
param_results1D = zeros(35,number_of_patients); % column i is the inferred parameters for patient i from minimizing cost over the active variables
param_results2D = zeros(35,number_of_patients);
costs1D = zeros(number_of_patients,1);
costs2D = zeros(number_of_patients,1);

for i = 1:2
    % patient data, cost function reference to each patient i
    patient_data = squeeze(patient_group_A(:,:,i));
    y_1 = w_active_1D'*random_parameter_array(1:35,i);
    y_2 = w_active_2D'*random_parameter_array(1:35,i);

    % functions to minimize, 1D and 2D active variables
    J_1 = @(x) cost2_eightCompModel(patient_data,IC,w_active_1D*x,HBs,pars.T);
    J_2 = @(x) cost2_eightCompModel(patient_data,IC,w_active_2D*x,HBs,pars.T);
    
    
    % results
    [xmin,fval,exitflag,output,lambda,grad,hessian] = fmincon(J_1,y_1,[],[],[],[],-1000,1000,[],options);
    [x2min,fval2,exitflag2,output2,lambda2,grad2,hessian2] = fmincon(J_2,y_2,[],[],[],[],[-1000,-1000],[1000,1000],[],options);
    
    results1D(i) = xmin;
    results2D(i,:) = x2min;
    
    param_results1D(:,i) = w_active_1D*xmin;
    param_results2D(:,i) = w_active_2D*x2min;
    
    costs1D(i) = fval;
    costs2D(i) = fval2;
end
    
save("optimizationResults","results1D","results2D","param_results1D","param_results2D");