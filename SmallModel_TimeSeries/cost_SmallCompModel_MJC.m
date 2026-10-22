function [y,qoi] = cost_SmallCompModel_MJC( IC,theta,HBs,T)

    % theta adult case for now
    theta(end+1) = T;
    [V,p,q,tplot] = GetModelOutput(IC,theta,HBs); % model output for new parameter set
    
    % vector for all the quantities of interest with the new parameter set
    qoi = zeros(52,1);
    
    % maximum and minimum volumes
    qoi(1) = max(V.LV);
    qoi(2) = min(V.LV);
    % maximum and minimum pressure in the systemic arteries
    qoi(3:52) = p.Ao; 
    
    y = sum(qoi);
    
end