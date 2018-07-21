function optz()
%     calc_optz('I:\study\Graduate\Summer\TraceAnalysis\1100-1149_readlink_sorted.txt')
%     calc_optz('I:\study\Graduate\Summer\TraceAnalysis\sorted_sample.txt')
    fid = fopen('I:\study\Graduate\Summer\TraceAnalysis\traces\MSNStorageCFS\IAT\new_filenames.txt');
    ln = fgetl(fid);
    mat=[];
    while ischar(ln)
        disp(ln);
%         data1 = load(ln)
        fnl_rsq = 0;
        fnl_jsd = 1;
        fnl_ks = 1;        
        fnl_para=[0 0 0 0];
        for n=1:5
            [para, resnorm, r_sq, ks_stat, jsd] = calc_optz(ln);        
            if r_sq > fnl_rsq
                %min_jsd = fval;                
                fnl_jsd = jsd;
                fnl_ks = ks_stat;
                fnl_rsq = r_sq;
                fnl_para= para;
            end
        end
        mat=[mat;fnl_jsd fnl_ks fnl_rsq fnl_para(1) 1-fnl_para(1) fnl_para(2) fnl_para(3) ];
        ln = fgetl(fid);                       
    end
    fclose(fid);
    csvwrite('output_lsq.csv',mat);
    fclose('all');
%     [para, resnorm] = calc_optz('I:\study\Graduate\Summer\TraceAnalysis\traces\MSNStorageCFS\IAT\CFS.2008-03-10.01-06.trace.csv.csv_HardFault.iat');    
%     t = "CFS.2008-03-10.01-06.trace.csv.csv_HardFault.iat"
%     mat=[resnorm para(1) 1-para(1) para(2) para(3) ]
    
    %calc_jsd('I:\study\Graduate\Summer\TraceAnalysis\1100-1149_readlink_sorted.txt', 3)
    %calc_jsd('I:\study\Graduate\Summer\TraceAnalysis\sorted_sample.txt', 3)
end

function [para, resnorm, r_sq, ks_stat, jsd]=calc_optz(fp)
    data1 = load(fp);
    [unqs, cdf_vals, pdf_vals] = nrm_data(data1);   
    %data1
    %x =[0,0]
    %cdf_vals, unqs = calc_cdf(data1);    
    %cdf_vals
    x0 = [rand(1,1), randi(10000,1, 2)];            
    cdf_vals = cdf_vals';
    pdf_vals = pdf_vals';
    %size(unqs);    
    %size(cdf_vals);
    %calc_lsq(x0, unqs, cdf_vals);
    %size(cdf_vals)
    f = @(x)calc_lsq(x, unqs, cdf_vals);
    lb = [0, 0, 0];
    ub = [1, inf, inf];
    size(lb);
    size(ub);
    size(x0);
%     options = optimoptions('lsqnonlin','Display','iter');
    [para, resnorm, residual] = lsqnonlin(f, x0, lb, ub);
    para;
    resnorm    ;
    y = cdf_vals;
    cdf_th = calc_cdf(para, unqs, 2, 0);
    r_sq = 1 - (sum((y - cdf_th).^2)/sum((y - mean(y)).^2));
    %para=para';    
    
    ks_stat = KS_stat(cdf_vals, cdf_th);   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    
%     x_axs = unique(data1);
%     plot(x_axs, cdf_th);
%     hold on
%     plot(x_axs, cdf_vals);
%     legend('Calculated', 'Actual');
%     title('Comparing CDFs')
%     xlabel('Data points')
%     ylabel('cdf')
    jsd = pre_jsd_calc(para, unqs, pdf_vals);
    
%     plot(cdf_vals);
%     hold on
%     plot(cdf_th);
%     legend('Data cdf', 'cdf th');
%     title('Actual vs Calculated');
%     xlabel('Points');
%     ylabel('Probability');
end

function f = pre_jsd_calc(x, unqs, pdf_vals)
    delta = 10^(-9);    
    intv_n = unqs - delta;
    intv_n(intv_n<0) = 0;
    intv_p = unqs + delta;    
%     intv_p(1) = 0;
    size(unqs);
    size(pdf_vals);
    cdf_n = calc_cdf(x, intv_n, 2, 0);
    cdf_p = calc_cdf(x, intv_p, 2, 0);
    pdf_calc = cdf_p - cdf_n;
    f = calc_jsd_div(pdf_vals, pdf_calc); 
end

function val = calc_jsd_div(P, Q)         
    size(P);
    size(Q);
    sum(P);    
    Q = Q./sum(Q);
    M = P+Q;
    M=M./2;
    sum(M);
%     find(isinf())
    val = (nansum(P.*log(P./M)) + nansum(Q.*log(Q./M)));%+ log(2).*(1-nansum(Q)));
    val = val/(2*log(2));
    % + log(2)*(1-nansum(Q));       
    %nansum(Q);
%     val = sum(P.*log(P./M)) + sum(Q.*log(Q./M));
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




function lsq = calc_lsq(x, data1, cdf_vals)
    prob = [x(1), 1-x(1)];    
    lmbd = x(2:3);
    size(prob);
    size(lmbd);
    size(data1);    
    t1 = (data1*lmbd);
    size(t1);
    t1 = -1*t1;
    t2 = exp(t1);
    t3 = t2*prob';    
    t3 = 1-t3;    
    lsq = t3 - cdf_vals;    
end


function cdf_th = calc_cdf(x, data1, k, f_min)
    if f_min
        prob = x(1:k);
        lmbd = x(k+1:2*k);
    else
        prob = [x(1), 1-x(1)];    
        lmbd = x(2:3);
    end
    
    size(prob);
    size(lmbd);
    size(data1);    
    t1 = (data1*lmbd);
    size(t1);
    t1 = -1*t1;
    t2 = exp(t1);
    t3 = t2*prob';    
    cdf_th = 1-t3;        
end

function stat = KS_stat(cdf1, cdf2)
    diff = cdf1 - cdf2;
    stat = max(abs(diff));
end


% function calc_jsd(fp, k)
%     data1 = load(fp);
%     [unqs, cdf_vals, pdf_vals] = nrm_data(data1);   
%     
%     %%%%%%%%%%%%% generate k random data%%%%%%
%     k=3;
%     r = randi(10, 1, k);
%     r = r./sum(r);
%     x0 = [r, randi(10,1, k)];
%     x0=x0'
%     size(x0);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5555
%     cdf_vals = cdf_vals';
%     [f, g]=jsd_grad(x0, 2, unqs, pdf_vals);
%     
% %     lb_pb = zeros(1,k);
% %     lb_lmbd = -1.*inf(1,k);
% %     ub_pb = ones(1,k);
% %     ub_lmbd = inf(1,k);
% %     lb = [lb_pb lb_lmbd];
% %     ub = [ub_pb ub_lmbd];
% %     lb=lb';
% %     ub = ub';
% %     A=[];
% %     b=[];
% %     Aeq = [ones(1,k) zeros(1,k)];
% %     beq = 1;
% %     nonlcon = [];
% %     f = @(x)jsd_grad(x,k, unqs, pdf_vals);
% %     
% %     options = optimoptions('fmincon','SpecifyObjectiveGradient',true);
% %     
% %     [para,fval] = fmincon(f, x0, A, b, Aeq, beq, lb, ub, nonlcon, options)    
% %     para = para';
% %     cdf_th = calc_cdf(para, unqs,k, 1);
% %     ks_1 = KS_stat(cdf_vals, cdf_th)    
% %     x_axs = unique(data1);
% %     plot(x_axs, cdf_th);
% %     hold on
% %     plot(x_axs, cdf_vals);
% %     legend('Calculated', 'Actual');
%     
% end
