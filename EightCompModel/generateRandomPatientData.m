%% generating randomized patient data

VSD_rad = 0.0;
HBs = 1;
% age = 7; unnecessary in the adult case
sex = 1; %F
percentile = 50;
weight = 70;
perturb = 0.5;
number_of_patients = 10;

% multiarray with volumes, pressures, and flows from each patient
patient_group_A = zeros(24,500,number_of_patients);

% mean values for each parameter
[IC, pars] = ComputeParameters_adult(weight,sex);
pars_array = makearray(pars);

% data from patients with randomized parameters
% normal distribution, mean = nominal value, variance = 5%*mean
random_patient_parameters = zeros(34,10);

for i=1:number_of_patients
    new_pars_array = (1+unifrnd(-perturb,perturb))*pars_array(1:end-1);
    % new_IC = normrnd(1,1.05)*IC;
    random_patient_parameters(:,i) = log(new_pars_array);
    
    % the time cycle is constant
    pars_struct = makestruct(new_pars_array,pars.T);
    [V,p,q,tplot] = GetModelOutput(IC,pars_struct,HBs); %newIC

    V = makearray(V);
    p = makearray(p);
    q = makearray(q);

    patient_group_A(:,:,i) = [V; p; q];
end    

save("randomizedPatientData","patient_group_A","random_patient_parameters","number_of_patients");