function jsd_pdf_single()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    fid = fopen('I:\study\Graduate\Summer\TraceAnalysis\traces\MSNStorageCFS\IAT\file_names.txt');
    ln = fgetl(fid);
    mat=[];
    while ischar(ln)
        disp(ln);
%         data1 = load(ln)
        fnl_jsd = 1;
        fnl_ks = 1;
        fnl_rsq = 0;
        fnl_para=[0];
        for n=1:5
            [para, fval, ks_stat, r_sq] = calc_jsd(ln)    
            if fval < fnl_jsd
                %min_jsd = fval;                
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
    csvwrite('jsdpdf_tmp.csv',mat);
    fclose('all');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    [para, fval, ks_stat] =  calc_jsd('I:\study\Graduate\Summer\TraceAnalysis\1100-1149_readlink_sorted.txt',2)
%     [para, fval, ks_stat] = calc_jsd('I:\study\Graduate\Summer\TraceAnalysis\sorted_sample.txt', 2)
end

function [para, fval, ks_stat, r_sq] = calc_jsd(fp)
    data1 = load(fp);
    [unqs, cdf_vals, pdf_vals] = nrm_data(data1);   

    %%%%%%%%%%%%% generate k random data%%%%%%    
    %r = randi(10, 1, k);
    %r = r./sum(r);
    x0 = randi(500);
    %[r, randi(10,1, k)];
%     x0=[0.9807 1-0.9807 238.7560 8.0838]
    %x0=x0';
    %size(x0);
%     jsd_pdf_grad(x0, k, unqs, pdf_vals')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %lb_pb = zeros(1,k);
    %lb_lmbd = zeros(1,k);
    %ub_pb = ones(1,k);
    %ub_lmbd = inf(1,k);
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
    %para = para';    
    fval;
    
    cdf_th = calc_cdf(para, unqs);
    size(cdf_vals);
    size(cdf_th) ;
    ks_stat = KS_stat(cdf_vals', cdf_th); 
    y = cdf_vals';
    r_sq = 1 - (sum((y - cdf_th).^2)/sum((y - mean(y)).^2));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     plot(cdf_vals');
%     hold on
%     plot(cdf_th);
%     legend('Data cdf', 'cdf th');
%     title('Actual vs Calculated');
%     xlabel('Points');
%     ylabel('Probability');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function cdf_th = calc_cdf(x, data1)    
    lmbd = x;    
    data1(1:3);
    t1 = (data1.*lmbd);
    t1(1:3);
    t1 = -1*t1;
    t2 = exp(t1);        
    cdf_th = 1-t2;        
end

function stat = KS_stat(cdf1, cdf2)
    diff = cdf1 - cdf2;
    stat = max(abs(diff));
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