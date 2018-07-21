function f = jsd_cdf_grad(x, k, unqs, cdf_vals)
    prob = x(1:k);
    lmbd = x(k+1:2*k);
    prob=prob';
    lmbd = lmbd';
    size(prob);
    size(lmbd);
    size(unqs);
    
    t1 = (unqs*lmbd);
    size(t1);
    t1 = -1*t1;
    t2 = exp(t1);
    t3 = t2*prob';    
    cdf_th = 1-t3;
    
    P = 1 - cdf_vals;
    P(end) = 0;
    Q = 1 - cdf_th;
    M = (P+Q)./2;
    
    lgP = log(P);
    P(size(P) -5:size(P));
    lgP(size(lgP) -5:size(lgP));
    lgP(end)
    P(end);
    
    
    lgQ = log(Q);
    Q(size(Q) -5:size(Q));
    lgQ(size(lgQ) -5:size(lgQ));
    lgQ(end);
    Q(end);
    
    f = -1*nansum(M.*log(M)) + (nansum(P.*log(P))/2 + nansum(Q.*log(Q)))/2
        
    g=[];
    
    
    
end