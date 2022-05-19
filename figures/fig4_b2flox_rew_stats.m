%% figure: b2flox
% (1) run fig1_achda_rew
% lag_wt = lag; val_wt = val; % save variables for comparison to b2flox

%% (2) extract lag,val for b2flox
beh = b2flox(19:30); % only reward recordings
% lag_b2 = lag; val_b2 = val;

%% STATS: b2flox
fig = figure; fig.Position(3) = 1000;
subplot(1,2,1); hold on
y = 2;
plot(lag_wt(:,y).*1000, val_wt(:,y), '.k', 'MarkerSize', 20);
plot(lag_b2(:,y).*1000, val_b2(:,y), '.m', 'MarkerSize', 20);
errorbar(nanmean(lag_wt(:,y).*1000), nanmean(val_wt(:,y)), SEM(val_wt(:,y),1), SEM(val_wt(:,y),1), SEM(lag_wt(:,y).*1000,1), SEM(lag_wt(:,y).*1000,1), '.k', 'MarkerSize', 20);
errorbar(nanmean(lag_b2(:,y).*1000), nanmean(val_b2(:,y)), SEM(val_b2(:,y),1), SEM(val_b2(:,y),1), SEM(lag_b2(:,y).*1000,1), SEM(lag_b2(:,y).*1000,1), '.m', 'MarkerSize', 20);
ylabel('DA peak amp'); ylim([0 20]); yticks([0:5:20]);
xlabel('time to rew (s)'); xlim([0 500]);
% p = ranksum(lag_wt(:,y),lag_b2(:,y)); p(2) = ranksum(val_wt(:,y),val_b2(:,y));
% title(sprintf('ranksum: lag = %1.3f, val = %1.3f',p(1),p(2))); axis('square');
[~,p] = ttest2(lag_wt(:,y),lag_b2(:,y)); [~,p(2)] = ttest2(val_wt(:,y),val_b2(:,y));
title(sprintf('ttest2: lag = %1.3f, val = %1.3f',p(1),p(2))); axis('square');

subplot(1,2,2); hold on
y = 3;
plot(lag_wt(:,y).*1000, val_wt(:,y), '.k', 'MarkerSize', 20);
plot(lag_b2(:,y).*1000, val_b2(:,y), '.g', 'MarkerSize', 20);
errorbar(nanmean(lag_wt(:,y).*1000), nanmean(val_wt(:,y)), SEM(val_wt(:,y),1), SEM(val_wt(:,y),1), SEM(lag_wt(:,y).*1000,1), SEM(lag_wt(:,y).*1000,1), '.k', 'MarkerSize', 20);
errorbar(nanmean(lag_b2(:,y).*1000), nanmean(val_b2(:,y)), SEM(val_b2(:,y),1), SEM(val_b2(:,y),1), SEM(lag_b2(:,y).*1000,1), SEM(lag_b2(:,y).*1000,1), '.g', 'MarkerSize', 20);
ylabel('ACh trough amp'); ylim([-8 0]); yticks([-8:2:0]);
xlabel('time to rew (s)'); xlim([200 700]);
% p = ranksum(lag_wt(:,y),lag_b2(:,y)); p(2) = ranksum(val_wt(:,y),val_b2(:,y));
% title(sprintf('ranksum: lag = %1.3f, val = %1.3f',p(1),p(2))); axis('square');
[~,p] = ttest2(lag_wt(:,y),lag_b2(:,y)); [~,p(2)] = ttest2(val_wt(:,y),val_b2(:,y));
title(sprintf('ttest2: lag = %1.3f, val = %1.3f',p(1),p(2))); axis('square');