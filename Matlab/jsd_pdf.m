function jsd_pdf()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    %read file which stores filepath of traces%
    fid = fopen('I:\study\Graduate\Summer\TraceAnalysis\traces\MSNStorageCFS\IAT\file_names.txt');
    ln = fgetl(fid);
    mat=[];
    while ischar(ln)
        disp(ln);
        fnl_jsd = 1;
        fnl_ks = 1;
        fnl_rsq = 0;
        fnl_para=[0 0 0 0];
        for n=1:5
            [para, fval, ks_stat, r_sq] = calc_jsd(ln, 2) ;
            if fval < fnl_jsd & r_sq>0
                %min_jsd = fval;                
                fnl_jsd = fval;
                fnl_ks = ks_stat;
                fnl_rsq = r_sq;
                fnl_para= para;
            end
        end
%         fnl_para
%         fnl_jsd
%         fnl_ks
%         fnl_rsq
        mat=[mat;fnl_jsd, fnl_ks, fnl_rsq, fnl_para ];
        ln = fgetl(fid);                       
    end
    fclose(fid);
    csvwrite('output_jsd_pdf.csv',mat);
    fclose('all');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    [para, fval, ks_stat] =  calc_jsd('I:\study\Graduate\Summer\TraceAnalysis\1100-1149_readlink_sorted.txt',2)
%     [para, fval, ks_stat] = calc_jsd('I:\study\Graduate\Summer\TraceAnalysis\sorted_sample.txt', 2)
end

function [para, fval, ks_stat, r_sq] = calc_jsd(fp, k)
    data1 = load(fp);
    [unqs, cdf_vals, pdf_vals] = nrm_data(data1);   

    %%%%%%%%%%%%% generate k random data%%%%%%    
    r = randi(10, 1, k);
    r = r./sum(r);
    l = randi(10000,1, k)
    x0 = [r, l];
%     x0=[0.9807 1-0.9807 238.7560 8.0838]
    x0=x0'
    size(x0);
%     jsd_pdf_grad(x0, k, unqs, pdf_vals')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lb_pb = zeros(1,k);
    lb_lmbd = zeros(1,k);
    ub_pb = ones(1,k);
    ub_lmbd = inf(1,k);
    lb = [lb_pb lb_lmbd];
    ub = [ub_pb ub_lmbd];
    lb=lb';
    ub = ub';
    A=[];
    b=[];
    Aeq = [ones(1,k) zeros(1,k)];
    beq = 1;
    nonlcon = [];
    f = @(x)jsd_pdf_grad(x,k, unqs, pdf_vals');
    
   % options = optimoptions('fmincon','SpecifyObjectiveGradient',true);
    
    options = optimoptions('fmincon','Display','iter');
    [para, fval] = fmincon(f, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);    
    para = para';
    size(para);
    fval;
    
    cdf_th = calc_cdf(para, unqs,k, 1);
    size(cdf_vals);
    size(cdf_th) ;
    ks_stat = KS_stat(cdf_vals', cdf_th); 
    y = cdf_vals';
    r_sq = 1 - (sum((y - cdf_th).^2)/sum((y - mean(y)).^2));
    cdf_th;
    
%     x_min = min(data1);
%     x_max = max(data1);
%     diff = x_max-x_min;
%     dt = (data1 - x_min)./diff;  
    
    %cdf_th = calc_cdf(para, dt,k, 1);
    size(data1);
    size(cdf_th);
    size(cdf_vals);
    cdf_vals = cdf_vals';
    %cdfplot(data1)    
    %size(h)
%     plot(cdf_vals);
%     hold on
%     plot(cdf_th);
%     legend('Data cdf', 'cdf th');
%     title('Actual vs Calculated');
%     xlabel('Points');
%     ylabel('Probability');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
%     lmbd = lmbd';
%     prob = prob';
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

% function val = calc_jsd_div(P, Q)         
%     size(P);
%     size(Q);
%     M = P+Q;
%     M=M./2;
% %     find(isinf())
%     nansum(P.*log(P./M)) + nansum(Q.*log(Q./M))    
%     val=0;
% %     val = sum(P.*log(P./M)) + sum(Q.*log(Q./M));
% end