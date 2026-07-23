% Compute active subspace for 8 compartment model
% Author: Ethan Mahan
% 5/18/2026

%% reference data
% data based on a 50th percentile adult female, no VSD

VSD_rad = 0.0; % cm
HBs = 1;
age = 7; % days
sex = 1;
percentile = 50;

[IC,pars] = ComputeParameters_adult(70,VSD_rad,sex);

[V,p,q,plot] = GetModelOutput(IC,pars,HBs);

% converting the parameters to an array for convenience later
parscell = struct2cell(pars);
pars_array = cell2mat(parscell); 

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
m = 1000;
h = 1e-8;
percent_width = 0.75;

L = (1-percent_width)*pars_array(1:35,1);
U = (1+percent_width)*pars_array(1:35,1);

random_parameter_array = zeros(35,m);

for i = 1:m
    random_parameter_array(:,i) = unifrnd(L,U,[35,1]);
end


%% computing C
% empty monte carlo approximation
C = zeros(35,35);

tic
parfor i = 1:m
    grad = zeros(35,1);
    theta = pars;
    rand_pars = pars;
    % put random parameters into a struct to pass to the model 
    % (not sampling a random VSD width)
    rand_pars_vector = random_parameter_array(:,i);

    rand_Rpars = [rand_pars_vector(1); rand_pars_vector(2); rand_pars_vector(3); rand_pars_vector(4); rand_pars_vector(5); rand_pars_vector(6); rand_pars_vector(7); rand_pars_vector(8); 1e10];
    rand_Cpars = [rand_pars_vector(10); rand_pars_vector(11); rand_pars_vector(12); rand_pars_vector(13)];
    rand_LApars = [rand_pars_vector(14); rand_pars_vector(15); rand_pars_vector(16); rand_pars_vector(17); rand_pars_vector(18); rand_pars_vector(19)];
    rand_LVpars = [rand_pars_vector(20); rand_pars_vector(21); rand_pars_vector(22); rand_pars_vector(23); rand_pars_vector(24)];
    rand_RApars = [rand_pars_vector(25); rand_pars_vector(26); rand_pars_vector(27); rand_pars_vector(28); rand_pars_vector(29); rand_pars_vector(30)];
    rand_RVpars = [rand_pars_vector(31); rand_pars_vector(32); rand_pars_vector(33); rand_pars_vector(34); rand_pars_vector(35)];

    rand_pars.Rpars = rand_Rpars;
    rand_pars.Cpars = rand_Cpars;
    rand_pars.LApars = rand_LApars;
    rand_pars.LVpars = rand_LVpars;
    rand_pars.RApars = rand_RApars;
    rand_pars.RVpars = rand_RVpars;
    rand_pars.T = pars.T;

    % cost function value at parameter sample i
    ref_cost = cost_eightCompModel(data,IC,rand_pars,HBs);
    
    % finite difference in each direction for 35 parameters
    for j = 1:35
        
        % 9th parameter relates to VSD radius, which is not being perturbed now
        if j~=9
            
            perturbed_pars_array = rand_pars_vector*1;

            % perturbed_pars_array = pars_array*1; % unperturb parameter j-1
            perturbed_pars_array(j) = perturbed_pars_array(j) + h; % perturb parameter j

            % put the parameters back into structure form
            Rpars = [perturbed_pars_array(1); perturbed_pars_array(2); perturbed_pars_array(3); perturbed_pars_array(4); perturbed_pars_array(5); perturbed_pars_array(6); perturbed_pars_array(7); perturbed_pars_array(8); 1e10];
            Cpars = [perturbed_pars_array(10); perturbed_pars_array(11); perturbed_pars_array(12); perturbed_pars_array(13)];
            LApars = [perturbed_pars_array(14); perturbed_pars_array(15); perturbed_pars_array(16); perturbed_pars_array(17); perturbed_pars_array(18); perturbed_pars_array(19)];
            LVpars = [perturbed_pars_array(20); perturbed_pars_array(21); perturbed_pars_array(22); perturbed_pars_array(23); perturbed_pars_array(24)];
            RApars = [perturbed_pars_array(25); perturbed_pars_array(26); perturbed_pars_array(27); perturbed_pars_array(28); perturbed_pars_array(29); perturbed_pars_array(30)];
            RVpars = [perturbed_pars_array(31); perturbed_pars_array(32); perturbed_pars_array(33); perturbed_pars_array(34); perturbed_pars_array(35)];

            theta.Rpars = Rpars;
            theta.Cpars = Cpars;
            theta.LApars = LApars;
            theta.LVpars = LVpars;
            theta.RApars = RApars;
            theta.RVpars = RVpars;

            theta.T = pars.T;

            % update jth gradient component using forward finite difference
            grad(j) = (cost_eightCompModel(data,IC,theta,HBs) - ref_cost)/h;
        else
            % skip update for param 9
        end

    end

    % add outer product to monte carlo approximation
    C = C + grad*grad';
    
end
toc

% divide by the number of samples
C = 1/m*C;

% remove row 9 and column 9 which are all zeros
newC = [C(1:8,:); C(10:35,:)];
newC = [newC(:,1:8) newC(:,10:35)];

%% active subspace
% eigendecomposition of C, W eigenvectors, w eigenevalues
[W, w] = eig(newC);

save("active_subspace_eightComp","W","w","newC","random_parameter_array");

%% Plotting
