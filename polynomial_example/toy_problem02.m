%% setup includes time domain, reference data, confidence interval parameter
clear; clc;
close all;
p = 4; % number of parameters
m = 40; % number of data points
x_domain = linspace(0,2,m);
n = 200; % number of random samples to be used in computation of active subspace
errors = normrnd(0,1,[1 m]);
true_parameters = [0.5;0.5;0.5;0.5];
data = f(true_parameters,x_domain) + errors;
X = [ones(m,1), x_domain', x_domain.^2', x_domain'.^7];
confidence_interval_alpha = 0.05;
% cost defined with respect to reference data
J = @(a) cost_toy_problem(a,data,x_domain);

%% compute active subspaces
h = 1e-6; % stepsize for gradient of J
C = zeros(4,4);
C_1 = zeros(3,3);
C_2 = zeros(2,2);
random_vector_set = zeros(4,n); 
L = [-1;-1;-1;-1]; % bounds for parameter samples
U = [1;1;1;1];
for i=1:n
    random_vector_set(:,i) = unifrnd(L,U);
end
% compute gradient of the cost function at n different points in parameter
% space, compute C matrix with outer product of gradients
for i=1:n
    a = random_vector_set(:,i);
    grad_J = zeros(4,1);
    for j=1:4
        new_a = a;
        new_a(j) = a(j)+h;
        grad_J(j) = (J(new_a) - J(a))/h;
    end
    C = C + grad_J*grad_J';
    C_1 = C_1 + grad_J(1:3)*grad_J(1:3)';
    C_2 = C_2 + grad_J(1:2)*grad_J(1:2)';
end

C = 1/n*C;
C_1 = 1/n*C_1;
C_2 = 1/n*C_2;
[W, lambda] = eig(C);
[W_remove1, lambda_remove1] = eig(C_1);
[W_remove2, lambda_remove2] = eig(C_2);

%% optimization
k = 10; % number of initializations
results_normal = zeros(4,k); % regular optimization over all 4 parameters
costs_normal = zeros(k,1);
results_act_sub1d = zeros(4,k); % optimization over 1D active subspace
costs_act_sub1d = zeros(k,1);
results_modified_act_sub = zeros(4,k); % optimization over the most influential parameter, with a 1D active subspace for the noninfluential parameters
costs_modified_act_sub = zeros(k,1);
results_act_sub2d = zeros(4,k); % optimization over 2D active subspace
costs_act_sub2d = zeros(k,1);


W_1_1D = W(:,4); % 1D active subspace
W_2_1D = W(:,1:3); % inactive variables for 1D active subspace
W_1_2D = W(:,3:4); % 2D active subspace
W_2_2D = W(:,1:2); % inactive variables for 2D active subspace

% W_mod_4D = W_remove1;
% M_4D = [W_mod_4D, [0;0;0]; [0,0,0,1]];

M_2D = [W_remove1(:,3), [0;0;0];0, 1]; % optimizing with M_2D*z in the model gives z_1 as an active variable for a_1 a_2 and a_3, and z_2 is passed into the model as a_4
W_2_modified_2D = [W_remove1(:,1:2); [0, 0]]; % inactive variables for modified active subspace, zeros in the 4th row since the 4th parameter is ignored in the active subspace

J_act_1D = @(alpha) cost_toy_problem(W_1_1D*alpha + W_2_1D*[0.5,0.5,0.5]',data,x_domain);
J_act_2D = @(alpha) cost_toy_problem(W_1_2D*alpha + W_2_2D*[0.5,0.5]',data,x_domain);
J_act_modified = @(z) cost_toy_problem(M_2D*z + W_2_modified_2D*[0.5,0.5]',data,x_domain);

for i=1:k
    guess = unifrnd(L,U); % random guess in the parameter interval
    guess1 = guess(4);
    guess2 = [guess(3), guess(4)]';
    guess3 = guess2;
    [amin,Jval] = fmincon(J,guess,[],[],[],[],L,U);
    [amin2,Jval2] = fmincon(J_act_1D,guess1,[],[],[],[],0,1);
    [amin3,Jval3] = fmincon(J_act_modified,guess2,[],[],[],[],[0;0],[1;1]);
    [amin4,Jval4] = fmincon(J_act_2D,guess2,[],[],[],[],[0;0],[1;1]);
    results_normal(:,i) = amin;
    costs_normal(i) = Jval;
    results_act_sub1d(:,i) = W_1_1D*amin2 + W_2_1D*[0.5,0.5,0.5]';
    costs_act_sub1d(i) = Jval2;
    results_modified_act_sub(:,i) = M_2D*amin3 + W_2_modified_2D*[0.5,0.5]';
    costs_modified_act_sub(i) = Jval3;
    results_act_sub2d(:,i) = W_1_2D*amin4 + W_2_2D*[0.5,0.5]';
    costs_act_sub2d(i) = Jval4;
end

%% confidence intervals
% OLS estimate for linear model \theta = (X^TX)^{-1}X^TY where Y is the data
% variance estimate \sigma = 1/(n-p)*R^TR where R = f(\theta) - Y
% X = [1 X X^2 X^7] \in R^{m \times 4}

V = inv(X'*X);
% linear ols estimate and variance estimate
OLS_estimate = V*X'*data'; % transposing since data is a row vector, also this estimate is not even close to the optimized case or correct parameters
R = (f(OLS_estimate,x_domain) - data)';
variance_estimate = 1/(m-p)*(R'*R);

% parameter and output confidence interval
c = tinv(0.975,m-p);
CI_width_half_param = c*sqrt(variance_estimate*diag(V));
% CI_width_half_output = c*sqrt(variance_estimate*X*diag(V)*X');
% PI_width_half = c*sqrt(variance_estimate*X*diag(V)*X' + variance_estimate);

% confidence interval for estimated mean (there is a 95% chance the mean is
% here!)
true_parameter_CI = [OLS_estimate-CI_width_half_param,OLS_estimate+CI_width_half_param];

% confidence interval and prediction interval for model output
c_pred_half = zeros(40,1);
p_pred_half = zeros(40,1);
for i=1:m
    c_pred_half(i) = c*sqrt(variance_estimate*X(i,:)*(diag(V).*X(i,:)'));
    p_pred_half(i) = c*sqrt(variance_estimate*X(i,:)*(diag(V).*X(i,:)') + variance_estimate);
end

compare_costs = [costs_normal(end),costs_modified_act_sub(end),costs_act_sub2d(end),costs_act_sub1d(end)];

relative_error_modified = sum(abs(true_parameters - results_modified_act_sub)./true_parameters);
relative_error_1d = sum(abs(true_parameters - results_act_sub1d)./true_parameters);
relative_error_2d = sum(abs(true_parameters - results_act_sub2d)./true_parameters);
relative_error_normal = sum(abs(true_parameters - results_normal)./true_parameters);

compare_errors = [relative_error_normal(end),relative_error_modified(end),relative_error_2d(end),relative_error_1d(end)];

%% activity scores
% standard active subspace
activity_scores = zeros(p,1);
for i=1:4
    activity_scores(i,1) = sum(diag(lambda)'*W(i,:).^2');
end

%% plotting
% data points, curves with estimated parameter values

figure(1)
hold on
vars = ["a_1","a_2","a_3","a_4"];
bar(log(activity_scores))
xticks(1:4)
xticklabels(vars)
ylabel("activity score")



graph_lambda = sort(diag(lambda),"descend");

figure(2)
hold on
vars = ["\lambda_1","\lambda_2","\lambda_3","\lambda_4"];
bar(log(diag(lambda)))
xticks(1:4)
xticklabels(vars)
ylabel("log eigenvalue")

figure(3)
hold on
vars = ["normal","modified","2D","1D"];
bar(compare_costs)
xticks(1:4)
xticklabels(vars)
ylabel("minimum cost")

figure(4);
hold on
bar(compare_errors)
xticks(1:4)
xticklabels(vars)
ylabel("relative parameter error")

% figure(1)
% hold on
% plot(x_domain,data,'.',LineStyle='none');
% plot(x_domain,f(results_normal(:,1),x_domain));
% 
% figure(11)
% hold on
% plot(x_domain,data,'o',LineStyle='none');
% plot(x_domain,f(results_act_sub1d(:,1),x_domain));
% 
% figure(12)
% hold on
% plot(x_domain,data,'o',LineStyle='none');
% plot(x_domain,f(results_act_sub2d(:,1),x_domain));
% 
% figure(13)
% hold on
% plot(x_domain,data,'o',LineStyle='none');
% plot(x_domain,f(results_modified_act_sub(:,1),x_domain));
% 
% y = f(OLS_estimate,x_domain);
% y1 = y' - p_pred_half;
% y2 = y' + p_pred_half;
% 
% figure(14)
% hold on
% plot(x_domain,data,'o',LineStyle='none');
% plot(x_domain,X*OLS_estimate);
% plot(x_domain,f(true_parameters,x_domain),"Color","k","LineWidth",2);
% plot(x_domain,y1);
% plot(x_domain,y2);
% fill([x_domain';flipud(x_domain')], [y1;flipud(y2)], [0.5 0.5 0.5], 'EdgeColor','none','FaceAlpha',0.3);

%% saving
% filename pmr_ab for x domain [a,b]

save("pmr_02","compare_errors","compare_costs","lambda","results_normal","results_modified_act_sub", ...
    "results_act_sub1d","results_act_sub2d","data","x_domain","true_parameters");

