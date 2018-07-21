
function jsd_pdf()
    %read file which stores filepath of traces
    fid = fopen('file_names.txt');
    ln = fgetl(fid);
    mat=[];
    
    while ischar(ln)
        disp(ln);
        fnl_jsd = 1;
        fnl_ks = 1;
        fnl_rsq = 0;
        fnl_para=[0 0 0 0];
        %iterate to find the most minimum value where r_sq metric is
        %positive
        for n=1:5
            [para, fval, ks_stat, r_sq] = calc_jsd(ln, 2) ;
            if fval < fnl_jsd & r_sq>0                
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
    csvwrite('output_jsd_pdf.csv',mat);
    fclose('all');
end


% calculates the parameters of a hyper-exponential distribution using
% Jensen-Shannon divergence(jsd) as a objective function.
%input: 
%   fp - path to trace file
%   k - number of single exponential distributions in a hyper exponential 
%output:
% para = 1 x 2k array containing k probabilities and their k lambdas
% fval = value of jsd at local minimum
% ks_stat = value of Kolmogorov–Smirnov test statistic between the actual and
% calculated distribution
%r_sq = r^2 between the actual and calculated distribution
function [para, fval, ks_stat, r_sq] = calc_jsd(fp, k)
    data1 = load(fp);
    
    %normalize data
    [unqs, cdf_vals, pdf_vals] = nrm_data(data1);   
      
    %generate k random probabilities
    r = randi(10, 1, k);
    r = r./sum(r);
    
    %generate a k x 1 vector of random numbers: initial value of lambda 
    l = randi(10000,1, k)
    x0 = [r, l];
    x0=x0'   
    
    %specify bounds
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
    
    % function handle to a function that calculates JSD
    f = @(x)jsd_pdf_grad(x,k, unqs, pdf_vals');
    
    %if gradient is calculated
    % options = optimoptions('fmincon','SpecifyObjectiveGradient',true);
    
    options = optimoptions('fmincon','Display','iter');
    
    
    [para, fval] = fmincon(f, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);    
    
    %resultant parameters of the hyperexponential distribution
    para = para';    
    
    %calculate cdf using the resultant parameters
    cdf_th = calc_cdf(para, unqs,k);
    
    %calcultate the K-S statistic
    ks_stat = KS_stat(cdf_vals', cdf_th); 
    y = cdf_vals';
    r_sq = 1 - (sum((y - cdf_th).^2)/sum((y - mean(y)).^2));    
    
    cdf_vals = cdf_vals';
    
    %plot data
%     plot(cdf_vals);
%     hold on
%     plot(cdf_th);
%     legend('Data cdf', 'cdf th');
%     title('Actual vs Calculated');
%     xlabel('Points');
%     ylabel('Probability');

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
    total = size(dt, 1);    
    
    %find unique points
    unqs = unique(x);
    unqs = (unqs - x_min)./diff;   
  
    % calculate cdf and pdf
    v = unique(x);
    v = vertcat(v, v(end) + 1);
    [cnts, edges] = histcounts(x, v);
    
    vals = cnts./total;
    pdf_vals = vals;    
    cdf_vals = cumsum(vals);            
end

%calculate cdf of hyper exponential distribution (1-sum(p.*l.exp(-data.*l)))
%input:
%x - parameters of hyper exponential in a 1 x k array probabilities first
%then corresponding lambdas
%data1 - data points normalized
% k - number of exponential distributions in a hyper expo
function cdf_th = calc_cdf(x, data1, k)
    prob = x(1:k);
    lmbd = x(k+1:2*k);        
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