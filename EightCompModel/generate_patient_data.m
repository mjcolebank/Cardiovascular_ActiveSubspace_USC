%% generate patient data

%% patient 1 (70 kg female, reference patient in active subspace file)
VSD_rad = 0.0; % cm
HBs = 1;
age = 7; % days
sex = 1;
percentile = 50;

[IC_patient_1,pars_patient_1] = ComputeParameters_adult(70,VSD_rad,sex);

[V_patient_1,p_patient_1,q_patient_1,plot] = GetModelOutput(IC_patient_1,pars_patient_1,HBs);

% put the parameters in a matrix for convenience later
parscell = struct2cell(pars_patient_1);
pars_array = cell2mat(parscell); 

volumescell = struct2cell(V_patient_1);
pressurescell = struct2cell(p_patient_1);
flowscell = struct2cell(q_patient_1);

volumes = cell2mat(volumescell);
pressures = cell2mat(pressurescell);
flows = cell2mat(flowscell);

% 25 by 500 matrix containing the reference quantities of interest at 500 time points
data_patient_1 = [volumes; pressures; flows];

save("patient_1","data_patient_1","V_patient_1","p_patient_1","q_patient_1","plot","pars_patient_1");

%% patient 2 (90 kg male)
VSD_rad = 0.0; % cm
HBs = 1;
age = 7; % days
sex = 2;
percentile = 50;

[IC_patient_2,pars_patient_2] = ComputeParameters_adult(90,VSD_rad,sex);

[V_patient_2,p_patient_2,q_patient_2,plot] = GetModelOutput(IC_patient_2,pars_patient_2,HBs);

% put the parameters in a matrix for convenience later
parscell = struct2cell(pars_patient_2);
pars_array = cell2mat(parscell); 

volumescell = struct2cell(V_patient_2);
pressurescell = struct2cell(p_patient_2);
flowscell = struct2cell(q_patient_2);

volumes = cell2mat(volumescell);
pressures = cell2mat(pressurescell);
flows = cell2mat(flowscell);

% 25 by 500 matrix containing the reference quantities of interest at 500 time points
data_patient_2 = [volumes; pressures; flows];

save("patient_2","data_patient_2","V_patient_2","p_patient_2","q_patient_2","plot","pars_patient_2");