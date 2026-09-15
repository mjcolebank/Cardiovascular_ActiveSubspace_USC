function p = LinearElastance(V,t,pars,T,t_off)
%  p = LinearElastance(V,t,pars,T,t_off)
% 
%  Description: 
%   Linear elastance model based on volume, time and associated parameters.
%
%  Inputs:
%   V     -- volume, e.g., LA, LV, RA and RV.
%   t     -- time variable
%   pars  -- cell array containing: Emax (mmHg/mL), Emin (mmHg/mL), Vd (mL), 
%            T_peak (s), T_relax (s)
%   T     -- heart rate
%   t_off -- optional time offset
%
%  Output:
%   p - pressure
% 
%  Author:
%   Mitchel J. Colebank <mcoleban@uci.edu> (May 2022)
% 
%  Refactored:
%   Alfonso.Limon@choc.org (Jan 2024)
%  -----------------------------------------------------------------------

% Update time
t = t + t_off;
t = mod(t,T);

% Parse associated parameters 
%parsCell = num2cell(pars);
[Emax, Emin, Vd, T_peak, T_relax] = pars{:};
T1 = T_peak;
T2 = T_peak+T_relax;
% Generate peak and relax elastance
rho_peak  = pi * t  / T_peak;
rho_relax = pi * (t - T_peak) / T_relax;%(T2 - T1);

E_peak  = (Emax - Emin) * (1 - cos(rho_peak )) / 2 + Emin;
E_relax = (Emax - Emin) * (1 + cos(rho_relax)) / 2 + Emin;

% Conditional indexes
ind_peak  = t <= T1;
ind_relax = (t > T1) & (t <= T2);

% Conditional update
E = Emin * ones(size(t));           % default
E(ind_peak)  = E_peak(ind_peak);    % if t(i) <= T_peak
E(ind_relax) = E_relax(ind_relax);  % if t(i) <= T_relax

% Pressure (elastance x volume)
p = E.*(V-Vd);

