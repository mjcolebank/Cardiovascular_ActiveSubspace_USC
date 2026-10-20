function [IC,pars] = ComputePerturbedParameters_adult(weight,VSD_rad,sex,parameter_number,h)

% parameter number is the index of the parameter to be perturbed
% h is the stepsize of the perturbation

% [IC,param] = ComputeParameters(weight,VSD_rad,T)
%
%  Description:
%   Provides parameters for  hemodynamic model, which includes Ventricular
%   Septal Defect (VSD). It calculates the parameter values using assumed
%   literature values or data provided as input
%
%  Inputs:
%   weight  - the bodyweight (in Kg) of the patient (assumed pediatric)
%   VSD_rad - the radius of the VSD (cm)
%   T       - the cardiac cycle length (s)
%
%  Output:
%   IC    - the initial conditions of the model
%   param - the parameters of the model
%
%  Author:
%   Mitchel J. Colebank <mcoleban@uci.edu>
%
%  Refactored:
%   Alfonso.Limon@choc.org (Jan 2024)
%  -----------------------------------------------------------------------
% Define argument block for optional parameters and default values
% NOTE: should add age in future uses


% For total volume, we use doi: 10.2471/BLT.10.080010 and scale based on
% weight
% For standard dimensions: https://academic.oup.com/ehjcimaging/article-abstract/16/3/233/2400086
if sex==2 %male
    tot_Volume = 5000;
    T = 60/70;
else %female
    tot_Volume = 4500;
    T = 60/80;
end

vla_stress_rat = 1.0;
vlv_stress_rat = 1.0;
vsa_stress_rat = 0.27;
vsv_stress_rat = 0.25;%0.075; % see https://link.springer.com/article/10.1007/s40140-014-0080-7?form=MG0AV3
vra_stress_rat = 1.0;
vrv_stress_rat = 1.0;
vpa_stress_rat = 0.58;
vpv_stress_rat = 0.11;

% All in ml
Vla_init = 0.01*tot_Volume*vla_stress_rat; % Left atrium (mL)
Vlv_init = 0.02*tot_Volume*vlv_stress_rat; % Left ventricle (mL)
Vsa_init = 0.15*tot_Volume*vsa_stress_rat; % Systemic arteries (mL)
Vsv_init = 0.60*tot_Volume*vsv_stress_rat; % Systemic veins (mL)
Vra_init = 0.01*tot_Volume*vra_stress_rat; % Right atrium (mL)
Vrv_init = 0.02*tot_Volume*vrv_stress_rat; % Right ventricle (mL)
Vpa_init = 0.09*tot_Volume*vpa_stress_rat; % Pulmonary arteries (mL)
Vpv_init = 0.10*tot_Volume*vpv_stress_rat; % Pulmonary veins (mL)


% Initial pressure estimates (mmHg)
% Note: this can be made specific to the data if necessary
if sex==1
    psa_max  = 110;
    psa_min = 70;
    psa_mean = (psa_max+2.*psa_min)./3;
    ppa_max  = 17;
    ppa_mean = 12; % Mean pulmonary arterial pressure (higher in newborn) (mmHg)
    ppa_min = 5;
    psv_mean = 8;  % Mean systemic venous pressure (mmHg)
    ppv_mean = 4; % Mean pulmonary venous pressure (mmHg)
    plv_max  = 111; % Systolic LV pressure (lower in newborn) (mmHg)
    plv_min  = 5;  % Diastolic LV pressure (mmHg)
    prv_max  = 18; % Systolic RV pressure (higher in newborn) (mmHg)
    prv_min  = 3;  % Diastolic RV pressure (mmHg)
    pla_max  = 13;  % Systolic LA pressure (mmHg)
    pla_min  = 3;  % Diastolic LA pressure (mmHg)
    pra_max  = 6;  % Systolic LA pressure (mmHg)
    pra_min  = 2;  % Diastolic LA pressure (mmHg)

    Edv_LV = 100;
    Esv_LV = 30;
    Edv_RV = Edv_LV;%110;
    Esv_RV = Esv_LV;%40;

    max_LA = 40;
    min_LA = 20;
    max_RA = 45;
    min_RA = 15;
else
    psa_max  = 120;
    psa_mean = 93.33; % Mean systemic arterial pressure (lower in newborn) (mmHg)
    psa_min = 80;
    ppa_max  = 20;
    ppa_mean = 13; % Mean pulmonary arterial pressure (higher in newborn) (mmHg)
    ppa_min = 8;
    psv_mean = 10;  % Mean systemic venous pressure (mmHg)
    ppv_mean = 4; % Mean pulmonary venous pressure (mmHg)
    plv_max  = 121; % Systolic LV pressure (lower in newborn) (mmHg)
    plv_min  = 5;  % Diastolic LV pressure (mmHg)
    prv_max  = 21; % Systolic RV pressure (higher in newborn) (mmHg)
    prv_min  = 3;  % Diastolic RV pressure (mmHg)
    pla_max  = 16;  % Systolic LA pressure (mmHg)
    pla_min  = 3;  % Diastolic LA pressure (mmHg)
    pra_max  = 8;  % Systolic LA pressure (mmHg)
    pra_min  = 2;  % Diastolic LA pressure (mmHg)

    Edv_LV = 120;
    Esv_LV = 40;
    Edv_RV = Edv_LV;%130;
    Esv_RV = Esv_LV;%50;

    max_LA = 55;
    min_LA = 25;
    max_RA = 55;
    min_RA = 25;
end

SV = mean([Edv_RV-Esv_RV,Edv_LV-Esv_LV]);
HR = 1./T;  % heart rate, beats/s
CO = SV*HR; % cardiac output, ml/s
% Paramters
Rmval      = 0.5/CO;  % Mitral valve resistance
Raval      = 0.1/CO;  % Aortic valve resistance
Rtval      = 0.5/CO;  % Tricuspid valve resistance
Rpval      = 0.1/CO;  % Pulmonary valve resistance

% Todo: Add an age dependent resistance function (doi: 10.1115/1.4062779)
Rs         = (psa_mean-psv_mean)/CO; % Resistance of the systemic capillaries
Rvc        = (psv_mean-pra_min)/CO; % Resistance of the vena cava
Rp         = (ppa_mean-ppv_mean)/CO; % Resistance of the pulmonary capillaries
Rpv        = (ppv_mean-pla_min)/CO; % Resistance of the pulmonary veins


PH_factor = 1.0; % This scales the pulmonary pressures to increase resistance
prv_max = prv_max*PH_factor; ppa_mean = ppa_mean*PH_factor; ppv_mean = ppv_mean*PH_factor;

%%
% Assume radius is similar to Poiseuille, except use:
% R = 8 mu L / pi r^4 ~ 8 mu/ pi r^3
% NOTE: this is not an ideal law; the resistance shrinks quickly for larger
% VSD radii. This may be more data driven


mu = 2.4002e-04; %Blood density in mmHg * s
rho = 1.055; % 1.35951 g/cm s^2 = 1 mmHg
if isempty(VSD_rad) || VSD_rad==0
    Rvsd = 1e10; %No opening
else

    %% Orifice flow argument
    Cd = 1e0;
    Kvsd = Cd*pi.*(VSD_rad.^2)*sqrt(2.*1333.22/rho); % last term is conversion factor
    % note: rho = g/cm^3, p = mmHg = 1.35951 g/cm s^2, so
    % sqrt(p/rho) = sqrt(g/cm s^2 * cm^3 / g) = sqrt(cm^2 / s^2) = cm/s
    % and pi *r^2 = cm^2, so kvsd*sqrt(p) = cm^3/s = mL/s
    % From Burkhoff paper
    %         Kvsd = 9.4540e-02*VSD_rad.^2;
    %         Kvsd = 9.4540e-02*VSD_rad.^2;
    %           Kvsd = 1*(pi.*VSD_rad.^2).*sqrt(2.*rho)./sqrt(1-0.9^4);
    % We want resistance for consistency
    Rvsd = 1./Kvsd;
end

Rpars = [Rmval;Raval;Rtval;Rpval;Rs;Rvc;Rp;Rpv;Rvsd];

% Csa = Vsa_init/psa_max;  % Systemic arterial compliance
Csv = Vsv_init/psv_mean;  % Systemic venous compliance
% Cpa = Vpa_init/ppa_max;  % Pulmonary arterial compliance
Cpv = Vpv_init/ppv_mean;  % Pulmonary venous compliance

Csa = SV/(psa_max-psa_min);  % Systemic arterial compliance
Cpa = SV/(ppa_max-ppa_min);  % Pulmonary arterial compliance


Cpars = [Csa;Csv;Cpa;Cpv];

% Max and Min elastance parameters, dead volume in the heart, and timing
% parameters describing the start, duration, and end of contraction (ask
% MJC for clarification)
Emax_LA    = pla_max/min_LA;%Vla_init;
Emin_LA    = (pla_max+pla_min)/(max_LA+min_LA);
Vd_LA      = 0;
T_peak_LA  = T*0.1;
T_relax_LA = T*0.10;
t_offset   = T*0.15;
LApars = [Emax_LA;Emin_LA;Vd_LA;T_peak_LA;T_relax_LA; t_offset];

Emax_RA    = pra_max/min_RA;
Emin_RA    = (pra_max+pra_min)/(max_RA+min_RA);
Vd_RA      = 0;%5;
T_peak_RA  = T*0.1;
T_relax_RA = T*0.1;
t_offset   = T*0.15;
RApars = [Emax_RA;Emin_RA;Vd_RA;T_peak_RA;T_relax_RA; t_offset];

Emax_LV    = plv_max/Esv_LV;
Emin_LV    = plv_min/Edv_LV;
Vd_LV      = 0;%5;
T_peak_LV  = T*0.25;
T_relax_LV = T*0.3;
LVpars = [Emax_LV;Emin_LV;Vd_LV;T_peak_LV;T_relax_LV];

Emax_RV    = prv_max/Esv_RV;
Emin_RV    = prv_min/Edv_RV;
Vd_RV      = 0;%5;
T_peak_RV  = T*0.25;
T_relax_RV = T*0.3;
RVpars = [Emax_RV;Emin_RV;Vd_RV;T_peak_RV;T_relax_RV];

% gather all the parameters into one vector

prams = [Rpars; Cpars; LApars; LVpars; RApars;RVpars];

%% add outside loop in driver over the parameter number

% perturb a particular parameter indexed by i (loop over i in driver)

prams(parameter_number) = prams(parameter_number) + h;


% put back all the values from the big parameter vector (now perturbed)
% into arrays with the same names as before

Rpars = [prams(1); prams(2); prams(3); prams(4); prams(5); prams(6); prams(7); prams(8); prams(9)];
Cpars = [prams(10); prams(11); prams(12); prams(13)];
LApars = [prams(14); prams(15); prams(16); prams(17); prams(18); prams(19)];
LVpars = [prams(20); prams(21); prams(22); prams(23); prams(24)];
RApars = [prams(25); prams(26); prams(27); prams(28); prams(29); prams(30)];
RVpars = [prams(31); prams(32); prams(33); prams(34); prams(35)];

pars.Rpars = Rpars;
pars.Cpars = Cpars;
pars.LApars = LApars;
pars.LVpars = LVpars;
pars.RApars = RApars;
pars.RVpars = RVpars;
pars.T = T;

IC =  [min_LA; Edv_LV; Vsa_init; Vsv_init;
    min_RA;Edv_RV; Vpa_init; Vpv_init];

end