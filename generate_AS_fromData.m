% This script will change the active subspaces using the stored gradients
% of the QoI's
% Reason: if J_i = (D_i - M_i(theta))/D_i, where i is the specific data
% point corresponding to a compartment output (e.g., valve flow), then
%   (d/d theta) J_i = - (1/D_i) dM_i / d theta
% So given data, D, load in the gradients of the mode output and scale them
% Since we typically compute J = sum ((D_i - M_i(theta))/D_i)^2, we can use
% existing gradients of M to our advantage

function [W,lambda] = generate_AS_fromData(data,cost_weights,param_ids)
load Samples_for_AS.mat grad_evals model_evals random_parameter_array
%%
grad_evals = grad_evals(:,param_ids,:);
random_parameter_array = random_parameter_array(param_ids,:);
%
[n_qoi,n_param,n_samp] = size(grad_evals);
qdata = zeros(14,1);
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
% maximum flow rates in the four valves mitral, aortic, tricuspid,
% pulmonary
qdata(9) = max(data(17,:));
qdata(10) = max(data(18,:));
qdata(11) = max(data(19,:));
qdata(12) = max(data(20,:));
% mean right atrial pressure
qdata(13) = mean(data(11,:));
% maximum pulmonary artery pressure
qdata(14) = max(data(15,:));

%% calculate residual
% Include a weight here in case you want to value one dataset over the
% other (e.g., valve flows might be noisy)
r = cost_weights.*(qdata - model_evals)./qdata;

%% Reconstruct covariance by looping over samples
C = zeros(n_param,n_param);

for j=1:n_samp
    thetaj = random_parameter_array(:,j)';
    Gj = squeeze(grad_evals(:,:,j));
    Glog = Gj.*thetaj;
    weights = -2.0.* (cost_weights.*r(:,j)./ qdata);
    grad_J = Glog'*weights;
    C = C + grad_J*grad_J';

    c_temp = grad_J*grad_J';
    corr_temp = c_temp.*0;
    for ii=1:n_param
        for iii=1:n_param
            corr_temp(ii,iii) = c_temp(ii,iii)./sqrt(c_temp(ii,ii)*c_temp(iii,iii));
        end
    end
end
C = C./n_samp;

%% Now sort eigenvalues and return active subspaces
[V,D] = eig(C);

lambda = diag(D);
[lambda,idx] = sort(lambda,'descend');
W = V(:,idx);


end