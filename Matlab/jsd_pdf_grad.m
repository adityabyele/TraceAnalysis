function f = jsd_pdf_grad(x, k, unqs, pdf_vals)
%     prob = x(1:k);
%     lmbd = x(k+1:2*k);
%     prob=prob';
%     lmbd = lmbd';
%     size(prob);
%     size(lmbd);
%     size(data1);

    size(unqs);
    size(pdf_vals);   
    %x = abs(min(unqs)); %in case of negative numbers
%     x=0;
%     for x1=sort(unqs)
%         x1
%         if(x1==0)
%             continue            
%         else
%             x=x1
%             break
%         end
%     end
    
%     n=0
%     while (floor(x*10^n)~=x*10^n)
%     n=n+1
%     end
%     n=n+1
    delta = 10^(-9);
    
    intv_n = unqs - delta;
    intv_n(intv_n<0) = 0;
    intv_p = unqs + delta;    
%     intv_p(1) = 0;
    cdf_n = calc_cdf(x, intv_n, 2, 1);
    cdf_p = calc_cdf(x, intv_p, 2, 1);
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
    lmbd = lmbd';
    prob = prob';
    t1 = (data1*lmbd);
    size(t1);
    t1 = -1*t1;
    t2 = exp(t1);
    t3 = t2*prob';    
    cdf_th = 1-t3;        
end
