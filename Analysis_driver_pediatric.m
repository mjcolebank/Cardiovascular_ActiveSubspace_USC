% 8 Compartment ODE model that includes ventricular septal defect
%
% Author: Mitchel J. Colebank
% Date: 5/24/2022
%%
clear; clc; %close all
%%
% The mathematical model consists of four cardiac compartments (atria and
% ventricles) and four arterial/venous compartments. The VSD is simulated
% through a resistance term that allows flow from the LV into the RV during
% systolic contraction and/or high LV diastolic pressure
VSD_rad = 0.0; % cm
HBs = 1;
age = 7; % days
sex = 1;
percentile = 50;

[IC,pars] = ComputeParameters_adult(70,VSD_rad,sex); %[IC, pars]

% This function then scales those parameters to what we believe would be
% their value in a pediatric case
[pars,IC] = new_allometric_scale(pars,IC,age,sex,percentile); %[pars, IC]

%%
tic
[V,p,q,tplot] = GetModelOutput(IC,pars,HBs);
toc

Volumes = struct2cell(V); % volume for 8 compartments at 500 times
Pressures = struct2cell(p); % pressure for 8 compartments at 500 times
flows = struct2cell(q); % flows accross 9 different resistors at 500 times

Volumes = cell2mat(Volumes);
Pressures = cell2mat(Pressures);
flows = cell2mat(flows);

reference_solution = [Volumes; Pressures; flows];

%% Sensitivites

h = 10^-8; % stepsize for local sensitivity analysis

sensitivities_multiarray = zeros(25,500,35); % declare storage for the sensitivities, sensitivities_multiarray(:,:,1) gives all the sensitivities with respect to parameter 1
cstep_multiarray = zeros(25,500,35); % storage for sensitivities with complex step

% adult case forward finite difference
for k=1:35
    [IC, pars] = ComputePerturbedParameters_adult_logscale(70,VSD_rad,sex,k,h); % perturbing log of parameter by h, then exponentiate exp(log(\theta) + h)
    % [pars,IC] = new_allometric_scale(pars,IC,age,sex,percentile);
    [newV,newp,newq,tplot] = GetModelOutput(IC,pars,HBs); % new solution for perturbed parameter set (slow)
    
    % want to take derivative with respect to the reference solution to get
    % a 25 by 500 matrix, put in sensitivities_multiarray(:,:,i)

    % convert model output to matrix (fast) perform matrix operation for
    % derivative (fast)

    newVolumes = struct2cell(newV); % volume for 8 compartments at 500 times
    newPressures = struct2cell(newp); % pressure for 8 compartments at 500 times
    newflows = struct2cell(newq); % flows accross 9 different resistors at 500 times

    newVolumes = cell2mat(newVolumes);
    Pressures = cell2mat(newPressures);
    newflows = cell2mat(newflows);

    new_solution = [newVolumes; newPressures; newflows];
    
    % forward finite difference
    sensitivities_multiarray(:,:,k) = (reference_solution - new_solution)./h;

    % old forward finite difference
    % sensitivities_V(i,:) = (newV.LV - ref_solution_V.LV)/h; % sensitivities_V(i) \approxeq \frac{\partialV}{\partial_{\theta_i}}
    % sensitivities_p(i,:) = (newp.LV - ref_solution_p.LV)/h; % its a structure so for now just do the LV portion we can look at the rest later
end

% adult case complex step method
for k=1:35
    [IC, pars] = ComputePerturbedParameters_adult(70,VSD_rad,sex,k,h*1i); % perturbing parameter by i*h
    % [pars,IC] = new_allometric_scale(pars,IC,age,sex,percentile);
    [newV,newp,newq,tplot] = GetModelOutput(IC,pars,HBs); % new solution for perturbed parameter set (slow)
    
    % want to take derivative with respect to the reference solution to get
    % a 25 by 500 matrix, put in sensitivities_multiarray(:,:,i)

    % convert model output to matrix (fast) perform matrix operation for
    % derivative (fast)

    newVolumes = struct2cell(newV); % volume for 8 compartments at 500 times
    newPressures = struct2cell(newp); % pressure for 8 compartments at 500 times
    newflows = struct2cell(newq); % flows accross 9 different resistors at 500 times

    newVolumes = cell2mat(newVolumes);
    Pressures = cell2mat(newPressures);
    newflows = cell2mat(newflows);

    new_solution = [newVolumes; newPressures; newflows];

    cstep_multiarray(:,:,k) = imag(new_solution)./h;
end


% infant aged 1 day through 7 days

for k = 1:7
    
end
%% Plotting
% Sensitivites
hold on
plot(tplot,sensitivities_V); % volume
plot(tplot,sensitivities_p); % pressure
% Pressure
figure(10*sex+percentile); clf; 

subplot(1,2,1); hold on;
plot(tplot,p.LA,'g','LineWidth',3);
plot(tplot,p.LV,'k','LineWidth',3);
plot(tplot,p.SA,'r','LineWidth',3);
plot(tplot,p.SV,'b','LineWidth',3);
ylabel('Pressure (mmHg)')
xlabel('Time (s)')
legend('LA','LV','Sys Art','Sys Vein')
grid on;
set(gca,'FontSize',20)

subplot(1,2,2); hold on;
plot(tplot,p.RA,'g','LineWidth',3);
plot(tplot,p.RV,'k','LineWidth',3);
plot(tplot,p.PA,'r','LineWidth',3);
plot(tplot,p.PV,'b','LineWidth',3);
ylabel('Pressure (mmHg)')
xlabel('Time (s)')
legend('RA','RV','Pulm Art','Pulm Vein')
grid on;
set(gca,'FontSize',20)

%% ventricular PV loop
figure(100*sex+percentile); 
subplot(1,2,1);hold on;
plot(V.LV,p.LV,'LineWidth',3);
subplot(1,2,2); hold on;
plot(V.RV,p.RV,'LineWidth',3);
ylabel('Pressure (mmHg)')
xlabel('Volume (mL)')
grid on;
set(gca,'FontSize',20);

%% atrial PV loop
figure(1000*sex+percentile); clf; hold on;
plot(V.LA,p.LA,'r','LineWidth',3);
plot(V.RA,p.RA,'b','LineWidth',3);
ylabel('Pressure (mmHg)')
xlabel('Volume (mL)')
grid on;
set(gca,'FontSize',20);




%% Flows
figure(10000*sex+percentile);clf; hold on;
plot(tplot,q.MVal,':r','LineWidth',3);
plot(tplot,q.AVal,'r','LineWidth',3);
plot(tplot,q.TVal,':b','LineWidth',3);
plot(tplot,q.PVal,'b','LineWidth',3);
plot(tplot,q.VSD,'-.k','LineWidth',3);
ylabel('Flow (mL/s)')
xlabel('Time (s)')
grid on;
set(gca,'FontSize',20);
legend('Mitral Valve','Aortic Valve','Tricuspid Valve','Pulmonic Valve','VSD');

disp('QP/QS ratio:')
disp(mean(q.PVal)./mean(q.AVal))
disp('Stroke Volume')
disp([max(V.LV)-min(V.LV),max(V.RV)-min(V.RV)])

VSD_A = pi.*(VSD_rad.^2);


% complex-step method
% log scale parameter perturbations $\frac{df}{d\log(\theta)}$
% store sensitivities in matrix (parameters \times (8 \times 500)) 
% (8 compartments)