function y = cost_eightCompModel_withSubsets(data,IC,theta_in,param_ids,param_remove,Wact,Winact,beta,HBs,T,cost_weights,plot_flag)
    % theta adult case for now
    % theta needs to be a vector for the optimizer
    % theta needs to be a structure for the model

    n_act      = size(Wact,2);
    n_inact    = size(Winact,2);
    n_normparam = length(param_remove);

    w_active_full = zeros(34,n_act); 
    w_inact_full = zeros(34,n_inact);
    w_active_full(param_ids,1:n_act) = Wact;
    w_inact_full(param_ids,1:n_inact) = Winact;

    theta_normal = theta_in(1:n_normparam);
    x_AS         = theta_in(n_normparam+1:end);

    theta_vec = w_active_full*x_AS+w_inact_full*beta;
    theta_vec(param_remove) = theta_normal;
    % Remove dead volume
    Vd_remove = [15 21 26 32];
    theta_vec(Vd_remove) = -inf;
    %EXPONENTIATE IF LOG_SCALE
    theta = makestruct(exp(theta_vec),10e8,T);

    [V,p,q,tplot] = GetModelOutput(IC,theta,HBs); % model output for new parameter set
    
    % vector for all the quantities of interest with the new parameter set
     qoi = zeros(14,1);
    qdata = zeros(14,1);
    
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
    qdata(6) = min(data(1,:));

    % maximum and minimum volume in the right ventricle
    qoi(7) = max(V.RV);
    qoi(8) = min(V.RV);
    qdata(7) = max(data(6,:));
    qdata(8) = min(data(6,:));

    % maximum and minimum volume in the right atrium
    % qoi(9) = max(V.RA);
    % qoi(10) = min(V.LA);
    % qdata(9) = max(data(5,:));
    % qdata(10) = min(data(5,:));

    % maximum flow rates in the four valves mitral, aortic, tricuspid,
    % pulmonary
    qoi(9) = max(q.MVal);
    qoi(10) = max(q.AVal);
    qoi(11) = max(q.TVal);
    qoi(12) = max(q.PVal);
    qdata(9) = max(data(17,:));
    qdata(10) = max(data(18,:));
    qdata(11) = max(data(19,:));
    qdata(12) = max(data(20,:));

    % mean right atrial pressure
    qoi(13) = mean(p.RA);
    qdata(13) = mean(data(11,:));

    % maximum pulmonary artery pressure
    qoi(14) = max(p.PA); 
    qdata(14) = max(data(15,:));

    % compute cost (2 norm of the difference of qois from new parameter set
    % and qois from data)

    % maybe should find a way to scale this

    y = sum((cost_weights.*(qdata-qoi)./(qdata)).^2);
    

    if plot_flag==1
        figure;
        subplot(2,2,1); hold on;
        plot(data(1,:),data(9,:),'k');
        plot(V.LA,p.LA,'r');
        xline(qdata(5:6),'--k');

        subplot(2,2,2); hold on;
        plot(data(2,:),data(10,:),'k');
        plot(V.LV,p.LV,'r');
        xline(qdata(3:4),'--k');
        

        subplot(2,2,3); hold on;
        plot(p.SA,'r');
        plot(p.RA,'c')
        plot(p.PA,'b');
        yline(qdata(1:2),'--r');
        yline(qdata(13),'--c');
        yline(qdata(14),'--b');

        subplot(2,2,4); hold on;
        plot(data(6,:),data(12,:),'k');
        plot(V.RV,p.RV,'b');
        xline(qdata(7:8),'--k');

        figure;
        subplot(2,2,1); hold on;
        plot(data(17,:),'k');
        plot(q.MVal,'r');
        yline(qdata(9),'--k')

        subplot(2,2,2); hold on;
        plot(data(18,:),'k');
        plot(q.AVal,'r');
        yline(qdata(10),'--k')

        subplot(2,2,3); hold on;
        plot(data(19,:),'k');
        plot(q.TVal,'r');
        yline(qdata(11),'--k')

        subplot(2,2,4); hold on;
        plot(data(20,:),'k');
        plot(q.PVal,'r');
        yline(qdata(12),'--k')


        
    end

end