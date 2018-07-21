function optz()
    %read file which stores filepath of traces
    fid = fopen('I:\study\Graduate\Summer\TraceAnalysis\traces\MSNStorageCFS\IAT\new_filenames.txt');
    ln = fgetl(fid);
    mat=[];
        
    while ischar(ln)
        disp(ln);
        fnl_rsq = 0;
        fnl_jsd = 1;
        fnl_ks = 1;        
        fnl_para=[0 0 0 0];
        %iterate to find the most max value of r_sq metric
        for n=1:5
            [para, resnorm, r_sq, ks_stat, jsd] = calc_optz(ln);        
            if r_sq > fnl_rsq                
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
end


% calculates the parameters of a hyper-exponential distribution using
% Least squares as a objective function.
%input: 
%   fp - path to trace file
%output:
% para = 1 x 3 array containing 1 probabilities and their 2 lambdas
% resnorm = norm of residual at local minimum
%r_sq = r^2 between the actual and calculated distribution
% ks_stat = value of Kolmogorov–Smirnov test statistic between the actual and
% calculated distribution
%jsd = value of jsd at local minimum
function [para, resnorm, r_sq, ks_stat, jsd]=calc_optz(fp)
    data1 = load(fp);
    %normalize data
    [unqs, cdf_vals, pdf_vals] = nrm_data(data1);       
    
    %randomly generate one probability and two lambdas
    %the second probability is calculated by doing 1 - this probability at
    %the end
    x0 = [rand(1,1), randi(10000,1, 2)];            
    cdf_vals = cdf_vals';
    pdf_vals = pdf_vals';    
    f = @(x)calc_lsq(x, unqs, cdf_vals);
    lb = [0, 0, 0];
    ub = [1, inf, inf];    
    
    %optimize for least squares as objective function
%     options = optimoptions('lsqnonlin','Display','iter');
    [para, resnorm, residual] = lsqnonlin(f, x0, lb, ub);    
    y = cdf_vals;
    
    %calc cdf of hyper expo
    cdf_th = calc_cdf(para, unqs);
    
    %calculate metrics
    r_sq = 1 - (sum((y - cdf_th).^2)/sum((y - mean(y)).^2));        
    ks_stat = KS_stat(cdf_vals, cdf_th);       
    jsd = pre_jsd_calc(para, unqs, pdf_vals);
    
    %plot cdf actual vs calculated
%     plot(cdf_vals);
%     hold on
%     plot(cdf_th);
%     legend('Data cdf', 'cdf th');
%     title('Actual vs Calculated');
%     xlabel('Points');
%     ylabel('Probability');
end


%return the jensen shannon divergence value between two distributions
%input:
%x - lambda
%unqs - data points normalized
%pdf_vals - pdf of data
function f = pre_jsd_calc(x, unqs, pdf_vals)
    delta = 10^(-9);    
    intv_n = unqs - delta;
    intv_n(intv_n<0) = 0;
    intv_p = unqs + delta;    
    cdf_n = calc_cdf(x, intv_n);
    cdf_p = calc_cdf(x, intv_p);
    pdf_calc = cdf_p - cdf_n;
    f = calc_jsd_div(pdf_vals, pdf_calc); 
end

%helper function to calculate jsd
function val = calc_jsd_div(P, Q)         
    size(P);
    size(Q);
    sum(P);    
    Q = Q./sum(Q);
    M = P+Q;
    M=M./2;
    sum(M);
    val = (nansum(P.*log(P./M)) + nansum(Q.*log(Q./M)));
    val = val/(2*log(2));    
end

%helper function to normalize data
%input:
% x - data
%output:
%unqs - unique data points
%cdf_vals - cdf of data
%pdf_vals - pdf of data 
function [unqs, cdf_vals, pdf_vals] = nrm_data(x)
    %normalize
    x_min = min(x);
    x_max = max(x);
    diff = x_max-x_min;
    dt = (x - x_min)./diff;   
    
    % calculate cdf and pdf
    total = size(dt, 1);    
    unqs = unique(x);
    unqs = (unqs - x_min)./diff;               
    v = unique(x);
    v = vertcat(v, v(end) + 1);
    [cnts, edges] = histcounts(x, v);    
    vals = cnts./total;
    pdf_vals = vals;    
    cdf_vals = cumsum(vals);                
end



%function to calculate least squares
% calculates the cdf of the distribtution and finds the difference
% with the cdf of data
% the output of this function is used by lsqnonlin which squares the
% differences and adds them
%input:
%x - 1x3 array of probability and 2 lambdas
%data1 - data points
%cdf_vals - cdf of data 
%output:
%difference between the cdf of data and cdf calculated
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


%calculate cdf of hyper exponential distribution (1-sum(p.*l.exp(-data.*l)))
%input:
%x - parameters of hyper exponential in a 1 x k array probabilities first
%then corresponding lambdas
%data1 - data points normalized
function cdf_th = calc_cdf(x, data1)
    prob = [x(1), 1-x(1)];    
    lmbd = x(2:3);    
    t1 = (data1*lmbd);    
    t1 = -1*t1;
    t2 = exp(t1);
    t3 = t2*prob';    
    cdf_th = 1-t3;        
end

%calculate KS statistic
function stat = KS_stat(cdf1, cdf2)
    diff = cdf1 - cdf2;
    stat = max(abs(diff));
end