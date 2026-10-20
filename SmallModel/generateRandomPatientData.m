%% generating randomized patient data

number_of_patients = 10;
perturb = 0.5;
HBs = 1;
% multiarray with volumes, pressures, and flows from each patient
patient_group_A = zeros(7,500,number_of_patients);

% mean values for each parameter
[IC, pars] = ComputeParameters_adult();
T = pars(end);
% data from patients with randomized parameters
% normal distribution, mean = nominal value, variance = 5%*mean
random_patient_parameters = zeros(12,number_of_patients);

for i=1:number_of_patients
    new_pars_array = (1+unifrnd(-perturb,perturb))*pars(1:end-1);
    % new_IC = normrnd(1,1.05)*IC;
    random_patient_parameters(:,i) = new_pars_array;
    pars_in = [new_pars_array,T];
    % the time cycle is constant
    [V,p,q,tplot] = GetModelOutput(IC,pars_in,HBs); %newIC

    V = makearray(V);
    p = makearray(p);
    q = makearray(q);

    patient_group_A(:,:,i) = [V; p; q];
end    

save("randomizedPatientData","patient_group_A","random_patient_parameters","number_of_patients");