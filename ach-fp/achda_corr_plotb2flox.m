fig = figure; fig.Position(3) = 1375;

subplot(1,3,1); hold on
shadederrbar(lags/Fs, nanmean(corr_wt,2), SEM(corr_wt,2), 'k');
shadederrbar(lags/Fs, nanmean(corr_b2,2), SEM(corr_b2,2), 'm');
legend({'control','','b2flox'})
xlabel('Lag (s)'); xlim([-1 1]);
ylabel('cross-correlation'); ylim([-0.8 0.4]); yticks([-0.8:0.2:1.2]);
title(sprintf('ACh/DA (wt = %d | b2 = %d)',size(corr_wt,2),size(corr_b2,2)));
axis('square');

[a,b] = min(corr_wt);
[c,d] = min(corr_b2);
subplot(1,3,2); hold on
plot_me = {a,c};
violinplot(plot_me);
errorbar(cellfun(@nanmean,plot_me),cellfun(@nanstd,plot_me)./sqrt(cellfun(@length,plot_me)),'.k','MarkerSize',20);
xticklabels({'control','b2flox'});
ylabel('max. coefficient'); ylim([-1 0]); yticks([-1:0.2:0]);
[~,p] = ttest2(plot_me{1},plot_me{2});
title(sprintf('ranksum: %1.4f, ttest2: %1.4f',ranksum(plot_me{1},plot_me{2}),p));
axis('square');

subplot(1,3,3); hold on
plot_me = {lags(b)/Fs*1000,lags(d)/Fs*1000};
violinplot(plot_me);
errorbar(cellfun(@nanmean,plot_me),cellfun(@nanstd,plot_me)./sqrt(cellfun(@length,plot_me)),'.k','MarkerSize',20);
xticklabels({'control','b2flox'});
ylabel('max. latency (ms)'); ylim([-200 400]);
[~,p] = ttest2(plot_me{1},plot_me{2});
title(sprintf('ranksum: %1.4f, ttest2: %1.4f',ranksum(plot_me{1},plot_me{2}),p));
axis('square');
movegui(fig,'center');