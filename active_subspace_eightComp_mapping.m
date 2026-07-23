%% mapping the cost function over the active variables

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

%% setup
% nominal values mapped to active subspace 2d, 3d
alpha_2D = w_active_2D'*pars_base;
beta_2D = w_inact_2D'*pars_base;

alpha_3D = w_active_3D'*pars_base;
beta_3D = w_inact_3D'*pars_base;

% vsd not included in active subspace analysis, need to set vsd
% resistance for the model to make sense
vsd_offset = zeros(35,1);
vsd_offset(9) =  1e10;

% patient data, only 1 patient
patient_data = squeeze(patient_group_A(:,:,1));

% cost functions with different dimensions of active subspace
J_2 = @(x2) cost2_eightCompModel(patient_data,IC,w_active_2D*x2+w_inact_2D*beta_2D+vsd_offset,HBs,pars.T); % x2 is a 2d vector
J_3 = @(x) cost2_eightCompModel(patient_data,IC,w_active_3D*x+w_inact_3D*beta_3D+vsd_offset,HBs,pars.T);

% 3d active subspace
% +- x direction causes problems
% +y direction is fine, -y direction causes problems
% +z direction is fine, -z direction causes problems
m = 25;
pointlist = zeros(26,26);
for i=0:m
    for j=0:m
        pointlist(i+1,j+1) = J_3(alpha_3D + [0.05;0.05*i;0.05*j]);
    end
end
%% another mapping

m = 51;
xcoords = linspace(alpha_2D(1)-0.05,alpha_2D(1)+0.05,m);
ycoords = linspace(alpha_2D(2)-0.05,alpha_2D(2)+0.05,m);
surface_point_list = zeros(m,m);

for i=1:m
    for j=1:m
        theta = [xcoords(i);ycoords(j)];
        surface_point_list(i,j) = J_2(theta);
    end
end
%% mapping 2D region around alpha
m = 25;
points = zeros(m+1,m+1);
for i=0:m
    for j=0:m
        points(i+1,j+1) = J_2(alpha_2D + [0.05*i;0.05*j]);
    end
end

% in [-1,1] direction, furthest valid return was at alpha + [-0.04,0.04]
% major issues in the [-1,0] and [0,-1] directions
% in [1,-1] direction, furthest valid return was at alpha + [0.03,-0.03]

num = 50;
morepoints = zeros(num+1,num+1);
for i = 0:num
    for j = 0:num
        morepoints(i+1,j+1) = J_2(alpha_2D + 0.01*[i*sqrt(3)/4 + j*0.25;j*sqrt(3)/4 - i*0.25]);
    end
end

%% list of points to map
m = 10;

% 2d points
x2_list = linspace(alpha_2D(1)*0.97,alpha_2D(1)*1.03,m);
y2_list = linspace(alpha_2D(2)*0.97,alpha_2D(2)*1.03,m);

% 3d points
% x3_list = linspace(alpha_2D(1)*0.95,alpha_2D(1)*1.05,m);
% y3_list = linspace(alpha_2D(2)*0.95,alpha_2D(2)*1.05,m);
% z3_list = linspace(alpha_3D(3)*0.95,alpha_3D(3)*1.05,m);

% values of cost function at the points
answers_J2 = zeros(m,m);
answers_J3 = zeros(m,m,m);

%% cost function evaluation, only 2d case for now
tic
for i = 1:m
    % for j=1:m
    %     answers_J2(i,j) = J_2([x2_list(i);y2_list(j)]);
    %     % for k=1:20
    %     %     answers_J3(i,j,k) = J_3([x3_list(i);y3_list(j);z3_list(k)]);
    %     % end
    % end
end
toc

save("active_subspace_map","answers_J2","answers_J3");