% FUNCTION LV_Ao.m
% This is an introductory function example for the course on cardiovascular
% physiology: theory and practice.
%
% Author: Mitchel J. Colebank
% Date:5/17/2022


function dY = SmallCompModel(t,y,param)
% Differential equations
Vlv = y(1);
Vao = y(2);

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
T         = param(13);

% Get LV and Ao pressure
plv = NonlinearElastance(Vlv,t,[Ees,Vlv_d,A,B,T_peak,T_relax,T]);
pao = Vao./Cao;

% Flows through valves and systemic arteries
qmv = max(((P_LeftAtr-plv)./Rmv),0);
qav = max(((plv-pao)./Rav),0);
qcap = (pao-P_SysCap)./Rsysart;

% Update differential equations
dVlv = qmv - qav;
dVsa = qav - qcap;
dY = [dVlv; dVsa];
if any(isnan(dY(:))) || Vlv==inf || Vao==inf
    disp('inf solution in ODE')
end
end