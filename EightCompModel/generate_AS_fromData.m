% This script will change the active subspaces using the stored gradients
% of the QoI's
% Reason: if J_i = (D_i - M_i(theta))/D_i, where i is the specific data
% point corresponding to a compartment output (e.g., valve flow), then
%   (d/d theta) J_i = - (1/D_i) dM_i / d theta
% So given data, D, load in the gradients of the mode output and scale them
% Since we typically compute J = sum ((D_i - M_i(theta))/D_i)^2, we can use
% existing gradients of M to our advantage

function [W,lambda] = generate_AS_fromData(data,cost_weights,param_ids)
% load Samples_for_AS_NewTiming.mat grad_evals model_evals random_parameter_array
load AS_8Comp_logparam_10_14_26.mat grad_evals model_evals random_parameter_array
%%
grad_evals = grad_evals(:,param_ids,:);
random_parameter_array = random_parameter_array(param_ids,:);
%
[n_qoi,n_param,n_samp] = size(grad_evals);
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

%% calculate residual
% Include a weight here in case you want to value one dataset over the
% other (e.g., valve flows might be noisy)
r = cost_weights.*(qdata - model_evals)./qdata;

%% Reconstruct covariance by looping over samples
C = zeros(n_param,n_param);

for j = 1:n_samp
    theta_j = random_parameter_array(:, j); % [n_param x 1]
    G_j     = squeeze(grad_evals(:, :, j)); % [n_qoi x n_param]
    M_j     = model_evals(:, j);           % [n_qoi x 1]
    
    % Residual vector: r_i = (qdata_i - M_i) / qdata_i
    % Cost J = sum_i [ cost_weights_i * r_i ]^2
    r_j = (qdata - M_j) ./ qdata;
    
    % Derivative dJ/dM_i = -2 * cost_weights_i^2 * (qdata_i - M_i) / (qdata_i^2)
    dJ_dM = -2.0 * (cost_weights.^2 .* r_j) ./ qdata; % [n_qoi x 1]
    
    % Chain rule: dJ/d_theta = G_j' * dJ_dM
    grad_theta_J = G_j' * dJ_dM;                      % [n_param x 1]
    
    % Chain rule for log-transformed parameters y = log(theta) => dJ/dy = theta * dJ/d_theta
    grad_log_J = theta_j .* grad_theta_J;             % [n_param x 1]
    
    % Accumulate unnormalized C matrix
    C = C + (grad_log_J * grad_log_J');
end
C = C./n_samp;

%% Now sort eigenvalues and return active subspaces
[V,D] = eig(C);

lambda = diag(D);
[lambda,idx] = sort(lambda,'descend');
W = V(:,idx);


end