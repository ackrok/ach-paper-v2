% mat_corr = struct;

% x = 1;
% mat_corr(x).fp1 = 'ACh-DLS';
% mat_corr(x).fp2 = 'DA-DLS';
% mat_corr(x).corr_imm = corr_an{1,1};

%%
for x = 1:length(mat_corr)
    [mat_corr(x).min_val, ii] = min(mat_corr(x).corr_imm([find(lags == 0):end],:));
    mat_corr(x).min_lag = lags(ii+find(lags == 0)-1)./50;
end

%%
fig = figure; fig.Position(3) = 1000; 
subplot(1,2,1); hold on; clr = {'g','b','r','m'};
plot([0 0],[-1 0.6],'--k');
for x = 1:length(mat_corr)
    plot(lags/Fs, nanmean(mat_corr(x).corr_imm,2), clr{x});
    %shadederrbar(lags/Fs, nanmean(mat_corr(x).corr_imm,2), SEM(mat_corr(x).corr_imm,2), clr{z});
end
shadederrbar(lags/Fs, nanmean(corr_an{1,3},2), nanmean(corr_an{1,2},2), 'k'); hold on
legend({'','ACh-DLS + DA-DLS','ACh-DMS + DA-DMS','ACh-DLS + DA-DMS','ACh-DMS + DA-DLS','',''});
xlabel('Lag (s)'); xlim([-1 1]);
ylabel('Correlation Coefficient'); ylim([-1 0.6]); yticks([-1:0.2:1]);
title(sprintf('ACh/DA corr IMM')); axis('square');

subplot(1,2,2); hold on;
a = reshape([mat_corr.min_val],[length(mat_corr(1).min_val) length(mat_corr)]);
% violinplot(a);
plot(a', 'k');
errorbar(nanmean(a,1), SEM(a,1),'.k', 'MarkerSize', 20);
xlim([0.5 4.5]); xticks([1:4]); xticklabels({'1','2','3','4'}); 
ylabel('Correlation Coefficient'); ylim([-1 0]); yticks([-1:0.2:1]);
% title(sprintf('coeff: signrank(i/m: %1.3f | i/r: %1.3f)', signrank(a(:,1),a(:,2)), signrank(a(:,1),a(:,3)))); axis('square');

% subplot(1,3,3); hold on;
% a = min_lag.*1000;
% violinplot(a);
% errorbar(nanmean(a,1), SEM(a,1),'.k', 'MarkerSize', 20);
% xlim([0.5 3.5]); xticks([1:3]); xticklabels({'imm','mov','rew'}); 
% ylabel('Latency to minimum (ms)'); ylim([0 250]); yticks([0:50:250]);
% title(sprintf('lag: signrank(i/m: %1.3f | i/r: %1.3f)', signrank(a(:,1),a(:,2)), signrank(a(:,1),a(:,3)))); axis('square');
% movegui(gcf,'center');

%%
b = nanmean(corr_an{1,2},1) - nanmean(corr_an{1,3},1);
b(4) = [];

p = [];
for x = 1:4
    [~,p(x)] = ttest(a(:,x),b(x));
end