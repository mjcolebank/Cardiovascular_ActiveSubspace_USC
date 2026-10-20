% FUNCTION NonlinearElastance.m
% This is an introductory function example for the course on cardiovascular
% physiology: theory and practice.
%
% Author: Mitchel J. Colebank
% Date:5/17/2022

function p = NonlinearElastance(V,t,pars)
%Input the following variables:
%t = time from 0 to T
%T = heart rate

Ees     = pars(1); % mmHg/mL
Vd      = pars(2); % mL
A       = pars(3); %1/ml
B       = pars(4); % mmHg
T_peak  = pars(5); % s
T_relax = pars(6); % s
T       = pars(7); % s


t   = mod(t,T);
% NOTE: this is similar in shape to elastance, but it is bounded in [0,1],
% hence we call it an 'activation function'
T1 = T_peak;
T2 = T_peak+T_relax;
% Generate peak and relax elastance
rho_peak  = pi * t  / T_peak;
rho_relax = pi * (t - T_peak) / T_relax;%(T2 - T1);

E_peak  = (1 - cos(rho_peak )) / 2;
E_relax = (1 + cos(rho_relax)) / 2;

% Conditional indexes
ind_peak  = t <= T1;
ind_relax = (t > T1) & (t <= T2);

% Conditional update
a = zeros(size(t));           % default
a(ind_peak)  = E_peak(ind_peak);    % if t(i) <= T_peak
a(ind_relax) = E_relax(ind_relax);  % if t(i) <= T_relax
% Pressure is a two part function; still dependent on volume
Pes = Ees.*(V-Vd);
Ped = B.*exp(A.*(V-Vd)) - B;
p = a.*(Pes-Ped) + Ped;

end

