function y = cost_SmallCompModel_withSubsets(data,IC,theta_in,param_ids,param_remove,Wact,Winact,beta,HBs,T,cost_weights,plot_flag)
    % theta adult case for now
    % theta needs to be a vector for the optimizer
    % theta needs to be a structure for the model
    n_par = 12;
    n_act      = size(Wact,2);
    n_inact    = size(Winact,2);
    n_normparam = length(param_remove);

    w_active_full = zeros(n_par,n_act); 
    w_inact_full = zeros(n_par,n_inact);
    w_active_full(param_ids,1:n_act) = Wact;
    w_inact_full(param_ids,1:n_inact) = Winact;

    theta_normal = theta_in(1:n_normparam);
    x_AS         = theta_in(n_normparam+1:end);

    theta_vec = w_active_full*x_AS+w_inact_full*beta;
    theta_vec(param_remove) = theta_normal;
    %EXPONENTIATE IF LOG_SCALE
    theta = [exp(theta_vec);T];

    [V,p,q,tplot] = GetModelOutput(IC,theta,HBs); % model output for new parameter set
    
    qdata = zeros(6,1);
    % maximum and minimum pressure in the systemic arteries
    qdata(1) = max(data(4,:));
    qdata(2) = min(data(4,:));
    % maximum and minimum volume in the left ventricle
    qdata(3) = max(data(1,:));
    qdata(4) = min(data(1,:));
    % maximum flow in the MV and AV valves
    qdata(5) = max(data(5,:));
    qdata(6) = max(data(6,:)); 

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

    y = sum((cost_weights.*(qdata-qoi)./(qdata)).^2);
    % y = sum((cost_weights.*(qdata-qoi)).^2);


    if plot_flag==1
        figure;
        subplot(1,3,1); hold on;
        plot(data(4,:),'k');
        plot(p.Ao,'r');
        yline(qdata(1:2),'--k');

        subplot(1,3,2); hold on;
        plot(data(1,:),'k');
        plot(V.LV,'r');
        yline(qdata(3:4),'--k');
        
        subplot(1,3,3); hold on;
        plot(data(5,:),'--b');
        plot(data(6,:),'--c');
        plot(q.MVal,'r');
        plot(q.AVal,'m');
        yline(qdata(5:6),'--k')


        
    end

end