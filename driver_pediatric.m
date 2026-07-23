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

% This function computes parameters and initial conditions for an adult
[IC,pars] = ComputeParameters_adult(70,VSD_rad,sex);

% This function then scales those parameters to what we believe would be
% their value in a pediatric case
[pars,IC] = new_allometric_scale(pars,IC,age,sex,percentile);

%%
tic
[V,p,q,tplot] = GetModelOutput(IC,pars,HBs);
toc

%% Plotting
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

