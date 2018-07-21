function jsd_cdf()
   calc_jsd('I:\study\Graduate\Summer\TraceAnalysis\1100-1149_readlink_sorted.txt',2)
%     [para, fval, ks_stat] = calc_jsd('I:\study\Graduate\Summer\TraceAnalysis\sorted_sample.txt', 2)
end

function calc_jsd(fp, k)
    data1 = load(fp);
    [unqs, cdf_vals, pdf_vals] = nrm_data(data1);   
    
    %%%%%%%%%%%%% generate k random data%%%%%%    
    r = randi(10, 1, k);
    r = r./sum(r);
    x0 = [r, randi(10,1, k)];
    x0=x0';
    size(x0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5555
    cdf_vals = cdf_vals';
    f =jsd_cdf_grad(x0, 2, unqs, cdf_vals);    
%     lb_pb = zeros(1,k);
%     lb_lmbd = -1.*inf(1,k);
%     ub_pb = ones(1,k);
%     ub_lmbd = inf(1,k);
%     lb = [lb_pb lb_lmbd];
%     ub = [ub_pb ub_lmbd];
%     lb=lb';
%     ub = ub';
%     A=[];
%     b=[];
%     Aeq = [ones(1,k) zeros(1,k)];
%     beq = 1;
%     nonlcon = [];
%     f = @(x)jsd_grad(x,k, unqs, pdf_vals);
%     
%     options = optimoptions('fmincon','SpecifyObjectiveGradient',true);
%     
%     [para,fval] = fmincon(f, x0, A, b, Aeq, beq, lb, ub, nonlcon, options)    
%     para = para';
%     cdf_th = calc_cdf(para, unqs,k, 1);
%     ks_1 = KS_stat(cdf_vals, cdf_th)    
%     x_axs = unique(data1);
%     plot(x_axs, cdf_th);
%     hold on
%     plot(x_axs, cdf_vals);
%     legend('Calculated', 'Actual');
    
end

function [unqs, cdf_vals, pdf_vals] = nrm_data(x)
    x_min = min(x);
    x_max = max(x);
    diff = x_max-x_min;
    dt = (x - x_min)./diff;   
    %%%%%%%%%%%%%%%%%%%%
    total = size(dt, 1);
    size(x);
    size(dt);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    unqs = unique(x);
    unqs = (unqs - x_min)./diff;   
    size(unqs);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %unqs = (unqs - x_min)./diff;   
    %size(unqs)
    v = unique(x);
    v = vertcat(v, v(end) + 1);
    [cnts, edges] = histcounts(x, v);
    %count=0
    %y = zeros(size(unqs))y = zeros(size(x));
    %for i = 1:length(x)
     %   y(i) = sum(x==x(i));
    %end
    %y
    %size(edges)
    %size(cnts)
    %edges(28875:size(edges))
    %unqs(28875:size(unqs))
    %cnts(28875:28876)
    %for i = 1:
     %   while 
    %size(cnts);
    vals = cnts./total;
    pdf_vals = vals;
    %cnts(1:5);
    cdf_vals = cumsum(vals);    
    %cdf_vals;
    
    size(cdf_vals);
end
