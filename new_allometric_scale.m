% allometrically scale parameters based on body weight (following from
% Hiebing et al.

function [pars,IC] = new_allometric_scale(pars,IC,age,sex,percentile)
b_BV = 1;
b_HR = -0.25;

%%
if sex==2 % Male
    if percentile == 25
        coeff = [1.2760e-08 -2.9576e-05 0.0264 3.3184];
    elseif percentile == 50
        coeff = [1.3426e-08 -3.1064e-05 0.0279 3.7063];
    elseif percentile == 75
        coeff = [1.4533e-08 -3.3328e-05 0.0297 4.0732];
    else 
        error('bad percentile');
    end
    tot_vol_0 = 5000; %(70Kg)

else % Female
    if percentile == 25
        coeff = [9.4973e-09 -2.3134e-05 0.0229 3.1651];
    elseif percentile == 50
        coeff = [1.0276e-08 -2.4571e-05 0.0242 3.5172];
    elseif percentile == 75
        coeff = [1.1178e-08 -2.6204e-05 0.0258 3.8549];
    else
        error('bad percentile');
    end
    tot_vol_0 = 4500; %(70Kg)
end


weight = polyval(coeff,age);
%%
% For total volume, we use doi: 10.2471/BLT.10.080010 and scale based on
% weight
weight0 = 70;
HR0 = 70;


tot_Volume = tot_vol_0.*(weight./weight0).^(b_BV); % ml/Kg for children over 1 month of age, else 100 ml/Kg?
HR = HR0.*(weight./weight0).^(b_HR);
pars.T = 60./HR;

% Distribute blood volume based on percentages
% Note: these percentages are the total STRESSED volume. Approximately 27%
% of the total blood volume is considered stressed in adults, but here we
% assume it is all stressed volume (need a reference later)
% All in ml
% Vla_init = 0.037*tot_Volume; % Left atrium (mL)
% Vlv_init = 0.088*tot_Volume; % Left ventricle (mL)
% Vsa_init = 0.110*tot_Volume; % Systemic arteries (mL)
% Vsv_init = 0.294*tot_Volume; % Systemic veins (mL)
% Vra_init = 0.037*tot_Volume; % Right atrium (mL)
% Vrv_init = 0.088*tot_Volume; % Right ventricle (mL)
% Vpa_init = 0.088*tot_Volume; % Pulmonary arteries (mL)
% Vpv_init = 0.258*tot_Volume; % Pulmonary veins (mL)

IC =  IC.*(weight./weight0).^(b_BV);

% MJC try to scale total volume
% vol_frac = IC./tot_vol_0.*1.1667;
% IC = vol_frac.*tot_Volume;
%% Apply allometric scaling
b_Ees = -1.0;
b_Emin = -0.5;%-0.5814;%0.7594;%-0.609;%-1.08;%-0.8;
b_V0 = 1;
b_tmax = -0.07;
b_R    =  -0.75;%-0.5;
b_valve = -0.25;
b_C = 1;
b_BV = 1;
b_pvc_svc = -1.333;

PVR0 = pars.Rpars(7);
SVR0 = pars.Rpars(5);
PVC0 = pars.Cpars(3);
SVC0 = pars.Cpars(1);

%%
% Rmv, Ra, Rt, Rp, Ras, Rvs, Rcp, Rap, Rvsd in VSD_circuit.png
valve_ids = [1:4];
% [Rmval, Raval, Rtval, Rpval, Rs, Rsv, Rp, Rpv, Rvsd] = pars.Rpars{:};
for j=valve_ids
    pars.Rpars(j) = pars.Rpars(j).*(weight./weight0).^(b_valve);
end

R_ids = [5:8];
% [Rmval, Raval, Rtval, Rpval, Rs, Rsv, Rp, Rpv, Rvsd] = pars.Rpars{:};
for j=R_ids
    pars.Rpars(j) = pars.Rpars(j).*(weight./weight0).^(b_R);
end
C_ids = [2 4];
for j=C_ids
    pars.Cpars(j) = pars.Cpars(j).*(weight./weight0).^(b_C);
end

%% PVR, PVC, SVR, SVC, and heart rate

c_pvr = 4;%7.1242;%10.6623;%4.76;%3.7;%2.91;%2.25;% %mmHg s/mL
c1_pvr = 1;%0.8167;%0.3748;%0.0642;%0.88;%1.3;%2.23;%2.62;%1.3;% days
c2_pvr = 2;%0.3139;%0.4312;%0.48;%1.6;%3.25;%4.29; days
Q_pvr = 1./(1.0+(age./c1_pvr).^c2_pvr);
% pars.Rpars(7) = c_pvr.*Q_pvr+pars.Rpars(7).*(1.0-Q_pvr);
pars.Rpars(7) = c_pvr.*Q_pvr+(PVR0.*(weight./weight0).^(-0.75)).*(1.0 - Q_pvr);
% pulm = c_pvr.*Q_pvr+.*(1.0-Q_pvr);
%% UPDATED SVR TO BE INVERSE
c_svr = 3.8;%12.4988;%4.7655;%2.92;%2.0;%5.0; %mmHg s/mL
c1_svr =50;%500;%183;%1000;%200%29.9;% days
c2_svr = 1.1;%0.7025;%0.7566;%0.12;%1.1;% days
Q_svr = 1./(1.0+(age./c1_svr).^c2_svr);
% pars.Rpars(5) = c_svr.*Q_svr+(SVR0.*(weight./weight0).^(-0.75)).*(1.0-Q_svr);


% c_svr = 2.5;%1.96;%2.06;%2.5;%0.5; %mmHg s/mL
% c1_svr = 400;%740;%484;%100;%90;% days
% c2_svr = 4.1;%0.12;%0.056;%0.25;%10.5;%4.1;% days
% Q_svr =  1./(1.0+(age./c1_svr).^c2_svr);
pars.Rpars(5) = c_svr.*Q_svr+(SVR0.*(weight./weight0).^(-0.5)).*(1.0-Q_svr);

% pars.Rpars(5) = (SVR0.*(weight./weight0).^(1.0));
% pars.Rpars(5) = c_svr.*Q_svr+pars.Rpars(5).*(1.0 - Q_svr);
% pars.Rpars(5) =c_svr + SVR0.*(weight./weight0).^(0.05);%.*(1-Q_svr);
% pars.Rpars(5) = c_svr.*Q_svr+(SVR0.*(weight./weight0).^(-0.25)).*(1.0 - Q_svr);
% pars.Rpars(5) = c_svr.*Q_svr+(SVR0.*(weight./weight0).^(-0.25));%.*(1.0 - Q_svr);

% sys = c_svr.*Q_svr+SVR0.*(1.0 - Q_svr);
%%
pars.Cpars(3) = PVC0.*(pars.Rpars(7)./PVR0).^b_pvc_svc;
pars.Cpars(1) = SVC0.*(pars.Rpars(5)./SVR0).^b_pvc_svc;

%% Cardiac parameters
% RApars = [Emax_RA;Emin_RA;Vd_RA;T_peak_RA;T_relax_RA; t_offset];
pars.RApars(1) = pars.RApars(1).*(weight./weight0).^(b_Ees);
pars.LApars(1) = pars.LApars(1).*(weight./weight0).^(b_Ees);
pars.RVpars(1) = pars.RVpars(1).*(weight./weight0).^(b_Ees);
pars.LVpars(1) = pars.LVpars(1).*(weight./weight0).^(b_Ees);

pars.RApars(2) = pars.RApars(2).*(weight./weight0).^(b_Emin);
pars.LApars(2) = pars.LApars(2).*(weight./weight0).^(b_Emin);
pars.RVpars(2) = pars.RVpars(2).*(weight./weight0).^(b_Emin);
pars.LVpars(2) = pars.LVpars(2).*(weight./weight0).^(b_Emin);

pars.RApars(4) = pars.T*0.10;%.*(weight./weight0).^(b_tmax);
pars.LApars(4) = pars.T*0.10;%.*(weight./weight0).^(b_tmax);
pars.RApars(5) = pars.T*0.20;%.*(weight./weight0).^(b_tmax);
pars.LApars(5) = pars.T*0.20;%.*(weight./weight0).^(b_tmax);
pars.RApars(6) = pars.T*0.15;%.*(weight./weight0).^(b_tmax);
pars.LApars(6) = pars.T*0.15;%.*(weight./weight0).^(b_tmax);

pars.RVpars(4) = pars.T*(0.25.*(weight./weight0).^(b_tmax));
pars.LVpars(4) = pars.T*(0.25.*(weight./weight0).^(b_tmax));
pars.RVpars(5) = pars.T*(0.55.*(weight./weight0).^(b_tmax));
pars.LVpars(5) = pars.T*(0.55.*(weight./weight0).^(b_tmax));



end