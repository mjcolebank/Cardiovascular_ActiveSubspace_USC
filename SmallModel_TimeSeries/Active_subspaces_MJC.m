% Compute active subspace for 8 compartment model
% Author: Ethan Mahan
% 5/18/2026

%% reference data
n_outputs = 52;
n_param = 12; % Exclude the cardiac cycle length
n_sample = 3000;
HBs = 1;
[IC,pars_OG] = ComputeParameters_adult();
T = pars_OG(end);
[V,p,q,~] = GetModelOutput(IC,pars_OG,HBs);

% converting the parameters to an array for convenience later
pars_array = pars_OG; 

volumescell = struct2cell(V);
pressurescell = struct2cell(p);
flowscell = struct2cell(q);

volumes = cell2mat(volumescell);
pressures = cell2mat(pressurescell);
flows = cell2mat(flowscell);

% 25 by 500 matrix containing the reference subject quantities of interest at 500 time points
data = [volumes; pressures; flows];
  
%% sampling parameters
% number of samples, step size for forward difference, % width of intervals
% for parameter values
h = 1e-6;
percent_width = 0.5;

L = (1-percent_width).*pars_array(1:n_param);
U = (1+percent_width).*pars_array(1:n_param);
% L(1) = 1; U(1) = 16;
% L(2) = 10; U(2) = 30;
% L(3) = 0.01; U(3) = 0.5;
% L(4) = 0.01; U(4) = 0.5;
% L(5) = 0.5; U(5) = 3.0;
% L(6) = 0.4; U(6) = 3.0;
% L(9) = 0.1; U(9) = 0.4;
% L(10) = 0.1; U(10) = 0.4;
% L(11) = 0.0; U(11) = 15;
% L(12) = 0.2; U(12) = 2.0;


random_parameter_array = zeros(n_param,n_sample);
cost_evals = zeros(n_sample,1);
model_evals = zeros(n_outputs,n_sample);
grad_evals  = zeros(n_outputs,n_param,n_sample);
for i = 1:n_sample
    random_parameter_array(:,i) = unifrnd(L,U,[1,12]);
end


%% computing C
% empty monte carlo approximation

tic
parfor i = 1:n_sample
    disp(i)
    grad = zeros(n_param,1);
    grad_qoi  = zeros(n_outputs,n_param);
    % put random parameters into a struct to pass to the model 
    % (not sampling a random VSD width)
    rand_pars_vector = random_parameter_array(:,i);

    % cost function value at parameter sample i
    [ref_cost,ref_qoi] = cost_SmallCompModel_MJC(IC,rand_pars_vector,HBs,T);
    cost_evals(i)    = ref_cost;
    model_evals(:,i) = ref_qoi;
    % cost_evals(i)    = log(ref_cost);
    % model_evals(:,i) = log(ref_qoi);
    
    % finite difference in each direction for 35 parameters
    for j = 1:n_param
        
       
            perturbed_pars_array = log(rand_pars_vector*1);
            
            % perturbed_pars_array = pars_array*1; % unperturb parameter j-1
            perturbed_pars_array(j) = perturbed_pars_array(j) + h; % perturb parameter j
            
            % update jth gradient component using forward finite difference
            [new_cost,new_qoi] = cost_SmallCompModel_MJC(IC,exp(perturbed_pars_array),HBs,T);
            
            grad(j) = (new_cost - ref_cost)/h;
            grad_qoi(:,j)  = (new_qoi-ref_qoi)/h;
            
            % grad(j) = (log(new_cost) - log(ref_cost))/h;
            % grad_qoi(:,j)  = (log(new_qoi)-log(ref_qoi))/h;
    end

    % add outer product to monte carlo approximation
    % if ~any(abs(grad_qoi(:))==inf) || ~any(isnan(grad_qoi(:)))
        grad_evals(:,:,i) = grad_qoi;
    % else
    %     2
    % end
end
toc

save("AS_AoTimeSeries_SmallModel_log_ODE23s_10_13_26","random_parameter_array","grad_evals",'cost_evals','model_evals');

%% Plotting
