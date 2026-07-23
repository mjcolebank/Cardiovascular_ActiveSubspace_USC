function [V,p,q,tplot] = GetModelOutput(IC,pars,HBs)


%%
% Starting and ending time
tstart = 0;
tend   = 30*pars.T; % 30 cycles
tspace = tstart:1e-4:tend;
options = odeset('RelTol',1e-8,'AbsTol',1e-8); % Necessary for the ODE solver

% Solve the model
% tic
sol = ode45(@EightCompModel,[tspace(1) tspace(end)],IC,options,pars);
% toc
% We solved the model for 30 heartbeats; we want to extract the last two
% for plotting
id = find(tspace>=(tend-HBs.*pars.T),1);
tplot = linspace(tspace(id),tspace(end),500); yout = deval(sol,tplot); yout = real(yout);

% Extract the solutions to the two differential equations
V.LA = yout(1,:); V.LV = yout(2,:);
V.SA = yout(3,:); V.SV = yout(4,:);
V.RA = yout(5,:); V.RV = yout(6,:);
V.PA = yout(7,:); V.PV = yout(8,:);

% Now, recompute the pressures and flows
p.LA = LinearElastance(V.LA, tplot, num2cell(pars.LApars(1:end-1)),pars.T,pars.LApars(end));
p.LV = LinearElastance(V.LV, tplot, num2cell(pars.LVpars),pars.T,0);
p.RA = LinearElastance(V.RA, tplot, num2cell(pars.RApars(1:end-1)),pars.T,pars.RApars(end));
p.RV = LinearElastance(V.RV, tplot, num2cell(pars.RVpars),pars.T,0);

p.SA = V.SA./pars.Cpars(1);
p.SV = V.SV./pars.Cpars(2);
p.PA = V.PA./pars.Cpars(3);
p.PV = V.PV./pars.Cpars(4);

% Use the 'max' operator to keep positive flows for valves
q.MVal  = max(((p.LA-p.LV)./pars.Rpars(1)),0);
q.AVal  = max(((p.LV-p.SA)./pars.Rpars(2)),0);
q.TVal  = max(((p.RA-p.RV)./pars.Rpars(3)),0);
q.PVal  = max(((p.RV-p.PA)./pars.Rpars(4)),0);

% Non valve flow
q.SYS = (p.SA-p.SV)./pars.Rpars(5);
q.PUL = (p.PA-p.PV)./pars.Rpars(6);
q.VC  = (p.SV-p.RA)./pars.Rpars(7);
q.PV  = (p.PV-p.LA)./pars.Rpars(8);
q.VSD  = sign(p.LV-p.RV).*sqrt(abs(p.LV-p.RV))./pars.Rpars(9);



tplot = tplot-tplot(1); % Set time to 0 at start
end