%return the jensen shannon divergence value between two distributions
%input:
%x - lambda
%unqs - data points normalized
%pdf_vals - pdf of data

function f = jsd_pdf_grad_single(x, unqs, pdf_vals)    

    %discretize exponential
    delta = 10^(-9);    
    intv_n = unqs - delta;
    intv_n(intv_n<0) = 0;
    intv_p = unqs + delta;     
    cdf_n = calc_cdf(x, intv_n);
    cdf_p = calc_cdf(x, intv_p);
    pdf_calc = cdf_p - cdf_n;
    
    
    f = calc_jsd_div(pdf_vals, pdf_calc);        
    
end


%calculate jsd between two pdfs
function val = calc_jsd_div(P, Q)         
    size(P);
    size(Q);
    sum(P);    
    Q = Q./sum(Q);
    M = P+Q;
    M=M./2;
    sum(M);
    val = (nansum(P.*log(P./M)) + nansum(Q.*log(Q./M)));%+ log(2).*(1-nansum(Q)));
    val = val/(2*log(2));    
end


%calculate cdf
function cdf_th = calc_cdf(x, data1)    
    lmbd = x;    
    data1(1:3);
    t1 = (data1.*lmbd);
    t1(1:3);
    t1 = -1*t1;
    t2 = exp(t1);        
    cdf_th = 1-t2;        
end
