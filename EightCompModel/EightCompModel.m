function dy = EightCompModel(t,y,pars)
% dy = EightCompModel(t,y,pars,Vsd_model)
% 
%  Description: 
%   Eight_Comp_Model simulates the dynamics of a circulatory system using
%   an electrical analog of the system with Ventricular Septal Defect (VSD).
%   It calculates the volumes in different compartments of the system based 
%   on the current time, state and associated parameters.
%
%  Inputs:
%   t    - time variable, used in the computation of pressures.
%   y    - state vector containing volumes in each compartment: 
%          LA, LV, SA, SV, RA, RV, PA, PV (refer to VSD_circuit.png)
%   pars - structure containing system parameters, e.g. resistances (Rpars),
%          compliances (Cpars), elastance parameters for different 
%          compartments (LApars, LVpars, RApars, RVpars), and period (T).
%   Vsd_model - type of pressure loss; 1=Poiseuille, 2=Bernoulli
%
%  Output:
%   dy   - The rate of change of volumes in each compartment.
%
%  Author:
%   Mitchel J. Colebank <mcoleban@uci.edu>
% 
%  Refactored:
%   Alfonso.Limon@choc.org (Jan 2024)
%  -----------------------------------------------------------------------

% % Check inputs
% arguments
%     t (1,1) double 
%     y (8,1) double %{mustBeNonnegative}
%     pars struct  % {adding field checks is slow see dev note below} 
%     Vsd_model double
% end

% Unpacking parameters
yCell  = num2cell(y);           
Rpars  = num2cell(pars.Rpars);  
Cpars  = num2cell(pars.Cpars);  
LApars = num2cell(pars.LApars);
LVpars = num2cell(pars.LVpars);
RApars = num2cell(pars.RApars);
RVpars = num2cell(pars.RVpars);
T      = pars.T;

% Rmv, Ra, Rt, Rp, Ras, Rvs, Rcp, Rap, Rvsd in VSD_circuit.png
[Rmval, Raval, Rtval, Rpval, Rs, Rsv, Rp, Rpv] = Rpars{:};

% LA, LV, SA, SV, RA, RV, PA, PV in VSD_circuit.png
[Vla, Vlv, Vsa, Vsv, Vra, Vrv, Vpa, Vpv] = yCell{:};

% Csa, Csv, Cpa, Cpv in VSD_circuit.png
[Csa, Csv, Cpa, Cpv] = Cpars{:};

% Calculate pressures in the arterial and venous components
psa = Vsa./Csa;
psv = Vsv./Csv;
ppa = Vpa./Cpa;
ppv = Vpv./Cpv;

% Compute pressures
pla = LinearElastance(Vla,t, LApars(1:end-1),T,LApars{end});
plv = LinearElastance(Vlv,t, LVpars,T,0);
pra = LinearElastance(Vra,t, RApars(1:end-1),T,RApars{end});
prv = LinearElastance(Vrv,t, RVpars,T,0);

% Flows through valves
qmval = max(((pla-plv)./Rmval),0);
qaval = max(((plv-psa)./Raval),0);
qtval = max(((pra-prv)./Rtval),0);
qpval = max(((prv-ppa)./Rpval),0);


% Non valve flow
qs  = (psa-psv)./Rs;
qp  = (ppa-ppv)./Rp;
qsv = (psv-pra)./Rsv;
qpv = (ppv-pla)./Rpv;

% Update differential equations
dVla = qpv   - qmval;
dVlv = qmval - qaval;
dVsa = qaval - qs;
dVsv = qs    - qsv;
dVra = qsv   - qtval; 
dVrv = qtval - qpval;
dVpa = qpval - qp;
dVpv = qp    - qpv;
dy   = [dVla; dVlv; dVsa; dVsv; dVra; dVrv; dVpa; dVpv];
