function temp
    gs = GlobalSearch;
    ms = MultiStart;
%     data1 = load('I:\study\Graduate\Summer\TraceAnalysis\1100-1149_readlink_sorted.txt');
% %     x_axs = unique(data1);
%     
% %     nrm_data = log(data1);
% %     x2_axs = unique(nrm_data);
% %     plot(data1)
% %     plot(nrm_data)
%     
%     avg = mean(data1)
%     std_dev = std(data1)
%     data = data1 - avg;
%     data = data./(std_dev*std_dev);
%     
%     total = size(data, 1);
%     v = unique(data);
%     v = vertcat(v, v(end) + 1);
%     [cnts, edges] = histcounts(data, v);
%     vals = cnts./total;
%     pdf_vals = vals;
%     unqs = unique(data);
    
%     [unqs, cdf_vals, pdf_vals] = nrm_data(data1);       
%     plot(unqs, pdf_vals)
%     plot(data)
%     data(1)
    
end


function [unqs, cdf_vals, pdf_vals] = nrm_data(x)
%     x_min = min(x);
%     x_max = max(x);
%     diff = x_max-x_min;
%     dt = (x - x_min)./diff;   
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
