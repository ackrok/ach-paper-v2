%% RUN achda_da2achpeak
good = [10:16,20,22:29,32:35,40:41,44:47]; % reward
beh = modAChDA(good);

beh = b2flox;

%% PLOT: STATS -- DA to ACh pause -- comparing peak magnitude, peak latency
da_b2 = da2peak_an{1}; ach_b2 = ach2peak_an{1};

da_wt = da2peak_an{1}; ach_wt = ach2peak_an{1};

%% MAXIMUM
[v_wt, ii] = max(da_wt(find(sta_time == 0):find(sta_time == 0.26),:)); % find MIN within range preceding ACh peak
t_wt = sta_time(ii + find(sta_time == 0) - 1);

[v_b2, ii] = max(da_b2(find(sta_time == 0):find(sta_time == 0.26),:)); % find MIN within range preceding ACh peak
t_b2 = sta_time(ii + find(sta_time == 0) - 1);

%% DIFF: MAX - MIN
[a, ii1] = min(da_wt(find(sta_time == -0.26):find(sta_time == 0),:)); % find MIN within range preceding ACh peak
clea[b, ii2] = max(da_wt(find(sta_time == 0):find(sta_time == 0.26),:)); % find MAX within range following ACh peak
v_wt = b - a;
t_wt = sta_time(ii2 + find(sta_time == 0) - 1);

[a, ii1] = min(da_b2(find(sta_time == -0.26):find(sta_time == 0),:)); % find MIN within range preceding ACh peak
[b, ii2] = max(da_b2(find(sta_time == 0):find(sta_time == 0.26),:)); % find MAX within range following ACh peak
v_b2 = b - a;
t_b2 = sta_time(ii2 + find(sta_time == 0) - 1);

%%
fig = figure; fig.Position([3 4]) = [1000 900];
subplot(2,2,1); hold on
a = {t_wt.*1000, t_b2.*1000};
violinplot(a); 
errorbar(cellfun(@nanmean,a),cellfun(@nanstd,a)./sqrt(cellfun(@length,a)),'.k','MarkerSize',20);
xticklabels({'control','b2flox'}); 
ylabel('latency to maximum (ms)'); ylim([0 300]); yticks([0:100:500]);
p = []; [~,p] = ttest2(a{1},a{2});
title(sprintf('DA to ACh peak - WT:%d b2:%d ms \n latency ranksum: %1.4f | ttest: %1.4f',round(nanmean(a{1})),round(nanmean(a{2})),ranksum(a{1},a{2}),p)); 
axis('square');

subplot(2,2,2); hold on
a = {v_wt,v_b2}; 
violinplot(a); 
errorbar(cellfun(@nanmean,a),cellfun(@nanstd,a)./sqrt(cellfun(@length,a)),'.k','MarkerSize',20);
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'WT','b2 flox'});
ylabel('difference amplitude (%dF/F)'); ylim([0 10]); yticks([0:2:10]);
axis('square');
p = []; [~,p] = ttest2(a{1},a{2});
title(sprintf('DA to ACh peak - WT:%1.2f b2:%1.2f dF/F \n max: ranksum: %1.4f | ttest: %1.4f',(nanmean(a{1})),(nanmean(a{2})),ranksum(a{1},a{2}),p)); 

movegui(gcf,'center');

% PLOT: AVERAGES for all groups
clr = {'k','r'};
% fig = figure; fig.Position([3 4]) = [1000 420];
subplot(2,2,3);
plot([0 0],[-5 3],'--k'); hold on; 
shadederrbar(sta_time, nanmean(da_wt,2), SEM(da_wt,2), clr{1}); hold on
shadederrbar(sta_time, nanmean(da_b2,2), SEM(da_b2,2), clr{2}); hold on
xlabel('Latency to ACh peak (s)'); ylabel('rDA1m (%dF/F)'); xlim([-1 1]);
axis('square')
title(sprintf('DA to ACh peak (n = %d mice)',nAn));

subplot(2,2,4);
plot([0 0],[-4 14],'--k'); hold on; 
shadederrbar(sta_time, nanmean(ach_wt,2), SEM(ach_wt,2), clr{1}); hold on
shadederrbar(sta_time, nanmean(ach_b2,2), SEM(ach_b2,2), clr{2}); hold on
xlabel('Latency to ACh peak (s)'); ylabel('ACh3.0 (%dF/F)'); xlim([-1 1]);
axis('square')
title(sprintf('ACh to ACh peak (n = %d mice)',nAn));
movegui(gcf,'center');