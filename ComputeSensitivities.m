% Compute sensitivities

clear; clc; % close all;

VSD_rad = 0.0; % cm
HBs = 1;
age = 7; % days
sex = 1;
percentile = 50;

[IC,pars] = ComputeParameters_adult(70,VSD_rad,sex); %[IC, pars]

[V,p,q,tplot] = GetModelOutput(IC,pars,HBs);

Volumes = struct2cell(V); % volume for 8 compartments at 500 times
Pressures = struct2cell(p); % pressure for 8 compartments at 500 times
flows = struct2cell(q); % flows accross 9 different resistors at 500 times

Volumes = cell2mat(Volumes); % (8 by 500)
Pressures = cell2mat(Pressures); % (8 by 500)
flows = cell2mat(flows); % (9 by 500)

reference_solution = [Volumes; Pressures; flows]; % (25 by 500), adult case

%% Sensitivites

h = 10^-8; % stepsize for local sensitivity analysis

sensitivities_multiarray = zeros(25,500,35); % declare storage for the sensitivities, sensitivities_multiarray(:,:,1) gives all the sensitivities with respect to parameter 1
cstep_multiarray = zeros(25,500,35); % storage for sensitivities with complex step

% tic
% % adult case forward finite difference
% for k=1:35
%     [IC, pars] = ComputePerturbedParameters_adult_logscale(70,VSD_rad,sex,k,h); % perturbing log of parameter by h, then exponentiate exp(log(\theta) + h)
%     % [pars,IC] = new_allometric_scale(pars,IC,age,sex,percentile);
%     [newV,newp,newq,tplot] = GetModelOutput(IC,pars,HBs); % new solution for perturbed parameter set (slow)
% 
%     % want to take derivative with respect to the reference solution to get
%     % a 25 by 500 matrix, put in sensitivities_multiarray(:,:,i)
% 
%     % convert model output to matrix (fast) perform matrix operation for
%     % derivative (fast)
% 
%     newVolumes = struct2cell(newV); % volume for 8 compartments at 500 times
%     newPressures = struct2cell(newp); % pressure for 8 compartments at 500 times
%     newflows = struct2cell(newq); % flows accross 9 different resistors at 500 times
% 
%     newVolumes = cell2mat(newVolumes);
%     newPressures = cell2mat(newPressures);
%     newflows = cell2mat(newflows);
% 
%     new_solution = [newVolumes; newPressures; newflows];
% 
%     % forward finite difference
%     sensitivities_multiarray(:,:,k) = (reference_solution - new_solution)./h;
% 
%     % old forward finite difference
%     % sensitivities_V(i,:) = (newV.LV - ref_solution_V.LV)/h; % sensitivities_V(i) \approxeq \frac{\partialV}{\partial_{\theta_i}}
%     % sensitivities_p(i,:) = (newp.LV - ref_solution_p.LV)/h; % its a structure so for now just do the LV portion we can look at the rest later
% end
% toc

% adult case complex step method
% h = 10^(-12) % smaller stepsize for complex step

% for k=1:35
%     [IC, pars] = ComputePerturbedParameters_adult(70,VSD_rad,sex,k,h*1i); % perturbing parameter by i*h
%     % [pars,IC] = new_allometric_scale(pars,IC,age,sex,percentile);
%     [newV,newp,newq,tplot] = GetModelOutput(IC,pars,HBs); % new solution for perturbed parameter set (slow)
% 
%     % want to take derivative with respect to the reference solution to get
%     % a 25 by 500 matrix, put in sensitivities_multiarray(:,:,i)
% 
%     % convert model output to matrix (fast) perform matrix operation for
%     % derivative (fast)
% 
%     newVolumes = struct2cell(newV); % volume for 8 compartments at 500 times
%     newPressures = struct2cell(newp); % pressure for 8 compartments at 500 times
%     newflows = struct2cell(newq); % flows accross 9 different resistors at 500 times
% 
%     newVolumes = cell2mat(newVolumes);
%     newPressures = cell2mat(newPressures);
%     newflows = cell2mat(newflows);
% 
%     new_solution = [newVolumes; newPressures; newflows];
% 
%     cstep_multiarray(:,:,k) = imag(new_solution)./h;
% end


% infant aged 1 day through 7 days
% need to perturb the parameters that are allometrically scaled for each
% day

[adult_IC,adult_pars] = ComputeParameters_adult(70,VSD_rad,sex);

ped_sensitivities = zeros(25,500,35,7); % 25 quantities of interest, 500 time points, 35 parameters, 7 ages

tic
for day = 1:7
    [pars,IC] = new_allometric_scale(adult_pars,adult_IC,day,sex,percentile);
    
    [V,p,q,tplot] = GetModelOutput(IC,pars,HBs); 

    Volumes = struct2cell(V); % volume for 8 compartments at 500 times
    Pressures = struct2cell(p); % pressure for 8 compartments at 500 times
    flows = struct2cell(q); % flows accross 9 different resistors at 500 times

    Volumes = cell2mat(Volumes); % (8 by 500)
    Pressures = cell2mat(Pressures); % (8 by 500)
    flows = cell2mat(flows); % (9 by 500)
    
    reference_solution_ped = [Volumes; Pressures; flows]; % reference solution for day day

    % perturb scaled parameters -> struct 2 cell and perturb in a loop, put
    % back into a structure manually (cut and paste from somewhere)
    % once we have perturbed parameters, get the sensitivities all at once
    % and put into a multi multiarray
    
    new_pars = struct2cell(pars);
    pars_array = cell2mat(new_pars);
    

    for k = 1:35
    
    perturbed_pars_array = pars_array*1; % unperturb parameter k-1
    perturbed_pars_array(k) = exp(log(pars_array(k)) + h); % perturb parameter k

    % put the parameters back into structure form
    Rpars = [perturbed_pars_array(1); perturbed_pars_array(2); perturbed_pars_array(3); perturbed_pars_array(4); perturbed_pars_array(5); perturbed_pars_array(6); perturbed_pars_array(7); perturbed_pars_array(8); perturbed_pars_array(9)];
    Cpars = [perturbed_pars_array(10); perturbed_pars_array(11); perturbed_pars_array(12); perturbed_pars_array(13)];
    LApars = [perturbed_pars_array(14); perturbed_pars_array(15); perturbed_pars_array(16); perturbed_pars_array(17); perturbed_pars_array(18); perturbed_pars_array(19)];
    LVpars = [perturbed_pars_array(20); perturbed_pars_array(21); perturbed_pars_array(22); perturbed_pars_array(23); perturbed_pars_array(24)];
    RApars = [perturbed_pars_array(25); perturbed_pars_array(26); perturbed_pars_array(27); perturbed_pars_array(28); perturbed_pars_array(29); perturbed_pars_array(30)];
    RVpars = [perturbed_pars_array(31); perturbed_pars_array(32); perturbed_pars_array(33); perturbed_pars_array(34); perturbed_pars_array(35)];

    pars.Rpars = Rpars;
    pars.Cpars = Cpars;
    pars.LApars = LApars;
    pars.LVpars = LVpars;
    pars.RApars = RApars;
    pars.RVpars = RVpars;
    
    pars.T = perturbed_pars_array(36);
    
    % pass structure to GetModelOutput
    
    [newV,newp,newq,tplot] = GetModelOutput(IC,pars,HBs);

    newVolumes = struct2cell(newV); % volume for 8 compartments at 500 times
    newPressures = struct2cell(newp); % pressure for 8 compartments at 500 times
    newflows = struct2cell(newq); % flows accross 9 different resistors at 500 times

    newVolumes = cell2mat(newVolumes);
    newPressures = cell2mat(newPressures);
    newflows = cell2mat(newflows);

    new_solution = [newVolumes; newPressures; newflows]; % perturbed solution for perturbed parameter k

    ped_sensitivities(:,:,k,day) = (new_solution - reference_solution_ped)./h; % sensitivity for all quantities of interest for parameter k with days old case
    end

end
toc

N = zeros(25,35,7); % norms for sensitivities of 25 QOI for 35 parameters on 7 different days
    for day = 1:7
        for i = 1:25
            vecplot = zeros(35,500);
            for n = 1:n
                vecplot(n) = A(i,:,k,day);
            end
            for k = 1:35
                sens_norm = sqrt(sum(vecplot.^2,2));
                N(i,:,day) = sens_norm;
            end
        end
    end

save('ped_sens_multiarray.mat','ped_sensitivities')
save('norms.mat','N')