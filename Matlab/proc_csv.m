function proc_csv()
res = csvread("output_lsq.csv");
res_single = csvread("output_lsq_single.csv");
obj = "LSq"
%nrm = res(:,1);

%nrm = sort(nrm)
jsd = sort(res(:,1));
ks = sort(res(:, 2));
r2 = sort(res(:, 3));
mean(jsd)
std(jsd)


jsd_single = sort(res_single(:,1));
ks_single = sort(res_single(:, 2));
r2_single = sort(res_single(:, 3));
mean(jsd_single)
std(jsd_single)

x_axs = [1:144];

plot(x_axs, r2);
hold on
plot(x_axs, r2_single);
legend('HyperExpo', 'Expo');
title(strcat(obj, ' as Objective'));
xlabel('Traces');
ylabel('r2:sorted');
%hold on
%plot(x_axs, r2);


% nrm_single = res_single(:, 1);
% sz = size(nrm_single);
% sz = sz(:,1);
% diff_ks = res_single(:,2)-res(:,2);
% count_neg = sum((diff_ks<0)==1)
% x_axs = 1:sz;
%diff = nrm_single - nrm;
% diff = nrm - nrm_single; %lsq/-
% diff;
% avg = mean(diff)
% std_dev = std(diff)
% mx = max(diff)
% mn = min(diff)
% count_neg = sum((diff<0)==1)
% count_zero = sum((diff==0)==1)
% intr_0_001 = sum((diff>0 & diff<0.001)==1)
% intr_001_01 = sum((diff>=0.001 & diff<0.01)==1)
% intr_01_1 = sum((diff>=0.01 & diff<0.1)==1)
% intr_1_2 = sum((diff>=0.1 & diff<0.2)==1)
% intr_2 = sum((diff>=0.2)==1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%pie%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %X = [intr_001_01 intr_01_1  count_neg  count_zero intr_0_001 intr_1_2 intr_2] %jsd/-
% %X = [count_neg intr_0_001 intr_001_01 intr_01_1 intr_1_2] %jsd lsq
% X = [intr_0_001 intr_001_01 intr_1_2 count_neg  intr_2] %lsq/-
% % labels={'<3rd decimal place', '3rd Decimal place', '1st decimal place' '>=0.2'}
% % legend(labels);
% h = pie(X)
% %title('JSD:Exp - HyperExp');%jsd/-
% %title('JSD:LSQ - JSD');%jsd lsq
% title('r^2:HyperExp - Exp');%;lsq/-
% hText = findobj(h, 'Type', 'text');
% percentValues = get(hText, 'String');
% %txt ={'>0.001 & <0.01: '; '>0.01 & <0.1: ';'<0: ';'=0: ';'>0 & <0.001: '; '>0.1 & <0.2: '; '>=0.2: '};%jsd/-
% txt ={'0> & <0.001: '; '>0.001 & <0.01: '; '>0.1 & <0.2: ';'<0: '; '>=0.2: '};%lsq/-
% %txt ={'<0: ';'<0.001: '; '>0.001 & <0.01: '; '>0.01 & <0.1: ';'>0.1 & <0.2: '};%jsd lsq
% combinedtxt = strcat(txt, percentValues);
% oldExtents_cell = get(hText, 'Extent');
% oldExtents = cell2mat(oldExtents_cell);
% hText(1).String = combinedtxt(1);
% hText(2).String = combinedtxt(2);
% hText(3).String = combinedtxt(3);
% hText(4).String = combinedtxt(4);
% hText(5).String = combinedtxt(5);
% %hText(6).String = combinedtxt(6);%jsd lsq = 5;jsd/-=7;lsq/-=5
% %hText(7).String = combinedtxt(7);
% newExtents_cell = get(hText,'Extent'); % cell array
% newExtents = cell2mat(newExtents_cell); % numeric array 
% width_change = newExtents(:,3)-oldExtents(:,3);
% signValues = sign(oldExtents(:,1));
% offset = signValues.*(width_change/2);
% textPositions_cell = get(hText,{'Position'}); % cell array
% textPositions = cell2mat(textPositions_cell); % numeric array
% textPositions(:,1) = textPositions(:,1) + offset; % add offset 
% 
% hText(1).Position = textPositions(1,:);
% hText(2).Position = textPositions(2,:);
% hText(3).Position = textPositions(3,:);
% hText(4).Position = textPositions(4,:);
% hText(5).Position = textPositions(5,:);
%hText(6).Position = textPositions(6,:);
%hText(7).Position = textPositions(7,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% size(res_single)
%     x_axs = unique(data1);
% bar(x_axs, diff);
% hold on
% scatter(x_axs, nrm_single);
% legend('HyperExpo', 'Expo');
% title('Difference in R^2');
% xlabel('File number');
% ylabel('Diff');
end