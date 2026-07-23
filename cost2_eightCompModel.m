function y = cost2_eightCompModel(data,IC,theta_vector,HBs,T)

    % theta adult case for now
    % theta needs to be a vector for the optimizer
    % theta needs to be a structure for the model

    theta = makestruct(theta_vector,10e8,T);

    [V,p,q,tplot] = GetModelOutput(IC,theta,HBs); % model output for new parameter set
    
    % vector for all the quantities of interest with the new parameter set
    qoi = zeros(16,1);
    qdata = zeros(16,1);
    
    % maximum and minimum pressure in the systemic arteries
    qoi(1) = max(p.SA); 
    qoi(2) = min(p.SA); 
    qdata(1) = max(data(13,:));
    qdata(2) = min(data(13,:));

    % maximum and minimum volume in the left ventricle
    qoi(3) = max(V.LV);
    qoi(4) = min(V.LV);
    qdata(3) = max(data(2,:));
    qdata(4) = min(data(2,:));

    % maximum and minimum volume in the left atrium
    qoi(5) = max(V.LA);
    qoi(6) = min(V.LA);
    qdata(5) = max(data(1,:));
    qdata(6) = min(data(2,:));

    % maximum and minimum volume in the right ventricle
    qoi(7) = max(V.RV);
    qoi(8) = min(V.RV);
    qdata(7) = max(data(6,:));
    qdata(8) = min(data(6,:));

    % maximum and minimum volume in the right atrium
    qoi(9) = max(V.RA);
    qoi(10) = min(V.LA);
    qdata(9) = max(data(5,:));
    qdata(10) = min(data(5,:));

    % maximum flow rates in the four valves mitral, aortic, tricuspid,
    % pulmonary
    qoi(11) = max(q.MVal);
    qoi(12) = max(q.AVal);
    qoi(13) = max(q.TVal);
    qoi(14) = max(q.PVal);
    qdata(11) = max(data(17,:));
    qdata(12) = max(data(18,:));
    qdata(13) = max(data(19,:));
    qdata(14) = max(data(20,:));

    % mean right atrial pressure
    qoi(15) = mean(p.RA);
    qdata(15) = mean(data(11,:));

    % maximum pulmonary artery pressure
    qoi(16) = max(p.PA); 
    qdata(16) = max(data(15,:));

    % compute cost (2 norm of the difference of qois from new parameter set
    % and qois from data)

    % maybe should find a way to scale this

    y = sum(((qoi - qdata)./(qdata)).^2);
    
    % figure()
    % hold on
    % plot(p.SA)
    % plot(data(13,:))

end