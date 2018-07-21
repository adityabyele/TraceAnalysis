function [f, g] = jsd_grad(x, k, data1, pdf_vals)
    prob = x(1:k);
    lmbd = x(k+1:2*k);
    prob=prob';
    lmbd = lmbd';
    size(prob);
    size(lmbd);
    size(data1);
    f=[];
    g=[];
    %%%%%%%%%%%%%%calc th cdf%%%%%%%%%%%%%%%%%%
    pdf_th = calc_pdf(data1, lmbd, prob);
%     t1 = (data1*lmbd);
%     size(t1);
%     t1 = -1*t1;
%     t2 = exp(t1);
%     t3 = t2*prob';    
%     cdf_th = 1-t3;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%calc JSD%%%%%%%%%%%%%%%%%%    
%     size(cdf_vals);
%     size(cdf_th);
%     t1 = log(cdf_vals);    
%     t1 = cdf_vals.*t1;    
%     log_qi = log(cdf_th);
%     log_qi(1);
%     cdf_th(1:4);
%     t2 = cdf_th.*log_qi;
%     p_q = (cdf_vals + cdf_th);
%     size(p_q);
%     t3 = (p_q).*log(p_q./2);
%     t4 = t1+t2-t3;        
%     f = nansum(t4);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     %%%%%%%%%%%%%%%%calc gradient%%%%%%%%%%%%%%
%     log_p_q = log(cdf_vals+cdf_th);    
%     trm1 = (log_qi - log_p_q + +log(2))./2;    
%     trm1(1);
%     %%%%%%%%%%%%%%%%%%calc_lmbd_grad%%%%%%%%%%%%%
%     t2 = (data1*lmbd);
%     size(t2);
%     t2 = -1*t2;
%     t3 = exp(t2);
%     t3 = t3.*prob;
%     t4 = data1.*t3;    
%     t5 = trm1.*t4;
%     t5 = nansum(t5, 1);    
%     size(t5);
%     grad_lmbd = t5;
%     
%     %%%%%%%%%%%%%%%%%%%%%calc_prob_grad%%%%%%%%%%
%     t2 = (data1*lmbd);
%     size(t2);
%     t2 = -1.*t2;
%     t3 = -1.*exp(t2); 
%     size(t3);
%     size(trm1);
%     t4 = trm1.*t3;    
%     size(t4);
%     r = isinf(t4);
%     find(r);
%     size(t4);
%     t5 = ~isinf(t4);    
%     t4 = t5.*t4;
%     size(t4);
%     t4 = nansum(t4,1);
%     size(t4);
%     %size(t4);
%     grad_prob = t4;
%     g_tmp = [grad_prob grad_lmbd]; 
%     g = g_tmp';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    
end

function pdf_th = calc_pdf(data1, lmbd, prob)
    t1 = (data1*lmbd);
    t1 = -1*t1;
    t2 = exp(t1);
    t3 = lmbd.*t2;
    t4 = prob.*t3;
    pdf_th = sum(t4, 2);
    size(pdf_th)
    sum(pdf_th)    
end
