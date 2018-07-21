% calculate jsd of a single exponential distribution
function jsd_pdf_single()
    %read file which stores filepath of traces%
    fid = fopen('file_names.txt');
    ln = fgetl(fid);
    mat=[];
    while ischar(ln)
        disp(ln);
        fnl_jsd = 1;
        fnl_ks = 1;
        fnl_rsq = 0;
        fnl_para=[0];
        %iterate to find the most minimum value where r_sq metric is
        %positive
        for n=1:5
            [para, fval, ks_stat, r_sq] = calc_jsd(ln)                
            if fval < fnl_jsd                
                fnl_jsd = fval;
                fnl_ks = ks_stat;
                fnl_rsq = r_sq;
                fnl_para= para;
           end
        end        
        mat=[mat;fnl_jsd, fnl_ks, fnl_rsq, fnl_para ];
        ln = fgetl(fid);                       
    end
    fclose(fid);
    
    %write all data to resultant csv
    csvwrite('jsdpdf_tmp.csv',mat);
    fclose('all');
end

% calculates the parameters of a exponential distribution using
% Jensen-Shannon divergence(jsd) as a objective function.
%input: 
%   fp - path to trace file
%output
% para = lambda
% fval = value of jsd at local minimum
% ks_stat = value of Kolmogorov–Smirnov test statistic between the actual and
% calculated distribution
%r_sq = r^2 between the actual and calculated distribution
function [para, fval, ks_stat, r_sq] = calc_jsd(fp)
    data1 = load(fp);
    [unqs, cdf_vals, pdf_vals] = nrm_data(data1);   

    %generaterandom initial value for lambda
    x0 = randi(500);    
    
    %bounding values
    lb = [0];
    ub = [inf];
    lb=lb';
    ub = ub';
    A=[];
    b=[];
    Aeq = [];
    beq = [];
    nonlcon = [];
    f = @(x)jsd_pdf_single_grad(x, unqs, pdf_vals');
    
   % options = optimoptions('fmincon','SpecifyObjectiveGradient',true);
    
    options = optimoptions('fmincon','Display','iter');
    [para, fval] = fmincon(f, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);               
    
    cdf_th = calc_cdf(para, unqs);    
    ks_stat = KS_stat(cdf_vals', cdf_th); 
    y = cdf_vals';
    r_sq = 1 - (sum((y - cdf_th).^2)/sum((y - mean(y)).^2));
    
    %plot data
%     plot(cdf_vals');
%     hold on
%     plot(cdf_th);
%     legend('Data cdf', 'cdf th');
%     title('Actual vs Calculated');
%     xlabel('Points');
%     ylabel('Probability');
end


%calculate cdf of exponential distribution
%input:
%x - lambda
%data1 - data points normalized
function cdf_th = calc_cdf(x, data1)    
    lmbd = x;    
    data1(1:3);
    t1 = (data1.*lmbd);
    t1(1:3);
    t1 = -1*t1;
    t2 = exp(t1);        
    cdf_th = 1-t2;        
end

%calculate KS statistic
function stat = KS_stat(cdf1, cdf2)
    diff = cdf1 - cdf2;
    stat = max(abs(diff));
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
    
    %find padf and cdf
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