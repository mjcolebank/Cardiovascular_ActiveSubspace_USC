function [IC,param] = ComputeParameters_adult()
% [IC,param] = ComputeParameters(weight,VSD_rad,T)
%
%  Description:
%   Provides parameters for  hemodynamic model, which includes Ventricular
%   Septal Defect (VSD). It calculates the parameter values using assumed
%   literature values or data provided as input
%

%  Author:
%   Mitchel J. Colebank <mcoleban@uci.edu>
%
% First, define some initial volumes for the differential equation
Vlv_init = 120; % LV Volume
Vao_init = 50;  % Aortic Volume

% Now define the parameters for the system
P_LA     = 6;    % left atrial pressure (mmHg)
P_SysCap = 20;   % systemic capillary pressure (mmHg)
Rmv      = 0.01; % resistance in the mitral valve (mmHg s/ml)
Rav      = 0.01; % resistance in the aortic valve (mmHg s/ml)
Rart     = 1.2;  % resistance of the systemic arteries/arterioles (mmHg s/ml)
Ees      = 1.8;  % End systolic elastance (mmHg/ml)
A        = 0.032;% Nonlinear elastance term (1/ml) 
B        = 0.10; % Nonlinear elastance term (mmHg/ml)
T_peak   = 0.25; % peak elastance (s)
T_relax  = 0.25;  % end of systole (s)
T        = 1.0;  % cardiac cycle length (s)
Vlv_d    = 10;   % LV dead volume (ml)
Cao      = 1.5;  % aortic compliance (ml s / mmHg)

% Stack all the parameters into a vector
param = [P_LA,P_SysCap,Rmv,Rav,Rart,Ees,A,B,T_peak,T_relax,Vlv_d,Cao,T];
IC = [Vlv_init,Vao_init];

end