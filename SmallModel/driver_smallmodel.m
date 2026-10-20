% This is an introductory example for the course on cardiovascular
% function: theory and practice.
%
% Author: Mitchel J. Colebank
% Date: 5/17/2022
%%
clear; clc; close all;
%%
% In this example, we will investigate a simple 'open loop' model of
% cardiovascular function
%%
% The mathematical model consists of a fixed, left atrial pressure source,
% a dynamics left ventricle (modeled by an elastance function), and a
% compliance system artery, with a constant, systemic capillary pressure
% source.

[IC,param] = ComputeParameters_adult();

% Starting and ending time
tstart = 0;
T      = param(end); % Length of a heartbeat
tend   = 30.*T; % 30 cycles

tspace = tstart:1e-4:tend;
options = odeset('RelTol',1e-6,'AbsTol',1e-8); % Necessary for the ODE solver


% Now solve our system of ODEs
[V,p,q,tplot] = GetModelOutput(IC,param,2);
P_LeftAtr = param(1);
P_SysCap  = param(2);
%% Plotting
% Pressure
figure(1); clf; hold on;
plot(tplot,tplot.*0+P_LeftAtr,'--k','LineWidth',3);
plot(tplot,p.LV,'k','LineWidth',3);
plot(tplot,p.Ao,'r','LineWidth',3);
plot(tplot,tplot.*0+P_SysCap,'--r','LineWidth',3);
ylabel('Pressure (mmHg)')
xlabel('Time (s)')
legend('LA','LV','Ao','Sys Art')
grid on;
set(gca,'FontSize',20)

%% PV loop
figure(2); hold on;
plot(V.LV,p.LV,'r','LineWidth',3);
ylabel('Pressure (mmHg)')
xlabel('Volume (mL)')
grid on;
set(gca,'FontSize',20);

%% Flows
figure(3); hold on;
plot(tplot,q.MVal,'r','LineWidth',3);
plot(tplot,q.AVal,'b','LineWidth',3);
plot(tplot,q.Sys,'k','LineWidth',3);
ylabel('Flow (mL/s)')
xlabel('Time (s)')
grid on;
set(gca,'FontSize',20);
legend('MV','AV','Sys Art')



