%return the jensen shannon divergence value between two distributions
%input:
%x - parameters of hyper exponential in a 1 x k array probabilities first
%then corresponding lambdas
% k - number of exponential distributions in a hyper expo
%unqs - data points normalized
%pdf_vals - pdf of data
function f = jsd_pdf_grad(x, k, unqs, pdf_vals)

    % discretize the hyperexponential distrbution
    delta = 10^(-9);    
    intv_n = unqs - delta;
    intv_n(intv_n<0) = 0;
    intv_p = unqs + delta;    
    
    %calculate cdf
    cdf_n = calc_cdf(x, intv_n, 2);
    cdf_p = calc_cdf(x, intv_p, 2);
    
    % calculate pdf of hyper exponential
    pdf_calc = cdf_p - cdf_n;
    
    %calculate jsd
    f = calc_jsd_div(pdf_vals, pdf_calc);        
    
end


%calculate jsd given pdfof 2 distributions
function val = calc_jsd_div(P, Q)         
    %normalize Q
    Q = Q./sum(Q);
    M = P+Q;
    M=M./2;    
    val = (nansum(P.*log(P./M)) + nansum(Q.*log(Q./M)));%+ log(2).*(1-nansum(Q)));
    val = val/(2*log(2));    
end



%calculate cdf
function cdf_th = calc_cdf(x, data1, k)    
    prob = x(1:k);
    lmbd = x(k+1:2*k);            
    lmbd = lmbd';
    prob = prob';
    t1 = (data1*lmbd);    
    t1 = -1*t1;
    t2 = exp(t1);
    t3 = t2*prob';    
    cdf_th = 1-t3;        
end
