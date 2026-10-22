% This script will change the active subspaces using the stored gradients
% of the QoI's
% Reason: if J_i = (D_i - M_i(theta))/D_i, where i is the specific data
% point corresponding to a compartment output (e.g., valve flow), then
%   (d/d theta) J_i = - (1/D_i) dM_i / d theta
% So given data, D, load in the gradients of the mode output and scale them
% Since we typically compute J = sum ((D_i - M_i(theta))/D_i)^2, we can use
% existing gradients of M to our advantage

function [W,lambda] = generate_AS_fromData_SVD(data,cost_weights,param_ids)
load AS_SmallModel_log_10_12_26.mat grad_evals model_evals random_parameter_array
%%
grad_evals = grad_evals(:,param_ids,:);
random_parameter_array = random_parameter_array(param_ids,:);
%
[n_qoi,n_param,n_samp] = size(grad_evals);
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

%% Use a "cutoff" for which parameter combinations to use based on the cost function
% threshold = 10;
% J = sum(((qdata - model_evals)./qdata).^2);
% [who,where] = sort(J,'ascend');
% ids = where(who>=threshold);
% J(ids) = [];
% %% Restrict residuals and gradients to samples below the cost threshold
% model_evals(:,ids) = [];
% grad_evals(:,:,ids) = [];
% random_parameter_array(:,ids) = [];
% n_samp = numel(where)-numel(ids);

%% calculate residual
% Include a weight here in case you want to value one dataset over the
% other (e.g., valve flows might be noisy)
r = cost_weights.*(qdata - model_evals)./qdata;

%% Reconstruct covariance by looping over samples
grad_storage = zeros(n_param,n_samp);

for j=1:n_samp
    thetaj = random_parameter_array(:,j)';
    Gj = squeeze(grad_evals(:,:,j));
    Glog = Gj;%.*thetaj;
    weights = -2.0.* (cost_weights.*r(:,j)./ qdata);
    grad_J = Glog'*weights;
    grad_storage(:,j) = grad_J;
end
grad_storage = grad_storage./sqrt(n_samp);

%% Now sort eigenvalues and return active subspaces
[U,S,V] = svd(grad_storage);

lambda = diag(S).^2;
[lambda,idx] = sort(lambda,'descend');
W = U(:,idx);


end