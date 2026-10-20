function [y,qoi] = cost_SmallCompModel_MJC( IC,theta,HBs,T)

    % theta adult case for now
    theta(end+1) = T;
    [V,p,q,tplot] = GetModelOutput(IC,theta,HBs); % model output for new parameter set
    
    % vector for all the quantities of interest with the new parameter set
    qoi = zeros(6,1);
    
    % maximum and minimum pressure in the systemic arteries
    qoi(1) = max(p.Ao); 
    qoi(2) = min(p.Ao); 
    
    % maximum and minimum volumes
    qoi(3) = max(V.LV);
    qoi(4) = min(V.LV);

    % maximum flow rates across the valves
    qoi(5) = max(q.MVal);
    qoi(6) = max(q.AVal);

    % compute cost (2 norm of the difference of qois from new parameter set
    % and qois from data)

    % maybe should find a way to scale this

    % y = sum(((qoi - qdata)./(qdata)).^2);

    % just do sum of model outputs for now
    
    y = sum(qoi);
    
end