function optz_single()
    fid = fopen('I:\study\Graduate\Summer\TraceAnalysis\traces\MSNStorageCFS\IAT\new_filenames.txt');
    ln = fgetl(fid);
    mat=[];
    while ischar(ln)
        disp(ln);
        fnl_rsq = 0;
        fnl_jsd = 1;
        fnl_ks = 1;        
        fnl_para=0;
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
%         data1 = load(ln)            
        mat=[mat;fnl_jsd fnl_ks fnl_rsq para];
        ln = fgetl(fid);                   
    end
    fclose(fid);
    csvwrite('output_lsq_single.csv',mat);
    fclose('all');
end

function [para, resnorm, r_sq, ks_stat, jsd] = calc_optz(fp)
    data1 = load(fp);
    [unqs, cdf_vals, pdf_vals] = nrm_data(data1);
    x0 = randi(10000);
    cdf_vals = cdf_vals';
%     calc_lsq(x0, unqs, cdf_vals);        
    f = @(x)calc_lsq(x, unqs, cdf_vals);
    lb = [0];
    ub = [inf];
%     para =[];
%     resnorm=[];
%     options = optimoptions('lsqnonlin','Display','iter');
    [para, resnorm, residual] = lsqnonlin(f, x0, lb, ub);
    para;
    resnorm;
    y = cdf_vals;
    r_sq = 1 - resnorm/sum((y - mean(y)).^2);
    
    
    lmbd = para;        
    t1 = (unqs.*lmbd);    
    t1 = -1*t1;
    t2 = exp(t1);        
    t3 = 1-t2;        
    ks_stat = KS_stat(cdf_vals, t3);
    pdf_vals = pdf_vals';
    jsd = pre_jsd_calc(lmbd, unqs, pdf_vals);    
end

function lsq = calc_lsq(x, data1, cdf_vals)    
    lmbd = x;    
    data1(1:3);
    t1 = (data1.*lmbd);
    t1(1:3);
    t1 = -1*t1;
    t2 = exp(t1);        
    t3 = 1-t2;    
    lsq = t3 - cdf_vals;    
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
    v = unique(x);
    v = vertcat(v, v(end) + 1);
    [cnts, edges] = histcounts(x, v);
    vals = cnts./total;
    pdf_vals = vals;
    cdf_vals = cumsum(vals);            
    size(cdf_vals);
end


function stat = KS_stat(cdf1, cdf2)
    diff = cdf1 - cdf2;
    stat = max(abs(diff));
end


function f = pre_jsd_calc(x, unqs, pdf_vals)
    %size(pdf_vals)
    delta = 10^(-9);    
    intv_n = unqs - delta;
    intv_n(intv_n<0) = 0;
    intv_p = unqs + delta;    
%     intv_p(1) = 0;
    size(unqs);
    size(pdf_vals);
    cdf_n = calc_cdf(x, intv_n);
    cdf_p = calc_cdf(x, intv_p);
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

function cdf_th = calc_cdf (para, unqs)
    lmbd = para;        
    t1 = (unqs.*lmbd);    
    t1 = -1*t1;
    t2 = exp(t1);        
    cdf_th = 1-t2; 
end