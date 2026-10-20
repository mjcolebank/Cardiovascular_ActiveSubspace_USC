function [V,p,q,tplot] = GetModelOutput(IC,param,HBs)

T         = param(13);

%%
% Starting and ending time
tstart = 0;
tend   = 30*T; % 30 cycles
tspace = tstart:1e-4:tend;
options = odeset('RelTol',1e-6,'AbsTol',1e-6); % Necessary for the ODE solver

% Solve the model
% tic
sol = ode23s(@SmallCompModel,[tspace(1) tspace(end)],IC,options,param);
% toc
% We solved the model for 30 heartbeats; we want to extract the last two
% for plotting
id = find(tspace>=(tend-HBs.*T),1);
tplot = linspace(tspace(id),tspace(end),500); yout = deval(sol,tplot); yout = real(yout);

%%
% Unpack parameters
P_LeftAtr = param(1);
P_SysCap  = param(2);
Rmv       = param(3);
Rav       = param(4);
Rsysart   = param(5);
Ees       = param(6);
A         = param(7);
B         = param(8);
T_peak    = param(9);
T_relax   = param(10);
Vlv_d     = param(11);
Cao       = param(12);
%%
% Extract the solutions to the two differential equations
V.LV = yout(1,:); V.Ao = yout(2,:);


% Get LV and Ao pressure
p.LV = NonlinearElastance(V.LV,tplot,[Ees,Vlv_d,A,B,T_peak,T_relax,T]);
p.Ao = V.Ao./Cao;

% Flows through valves and systemic arteries
q.MVal = max(((P_LeftAtr-p.LV)./Rmv),0);
q.AVal = max(((p.LV-p.Ao)./Rav),0);
q.Sys = (p.Ao-P_SysCap)./Rsysart;

tplot = tplot-tplot(1); % Set time to 0 at start
% figure(1);clf; plot(q.AVal); hold on; plot(q.MVal);
end