%% PLOT N = X mice
fig = figure; fig.Position(3) = 1375; 
lbl = {'immobility','locomotion','reward'}; clr = {'r','g','b'};
for z = 1:3
    subplot(1,3,z); hold on
    plot([0 0],[-1 0.5],'--k');
    a = nanmean(corr_an{z,3},2);
    a = a - nanmean(a(find(lags/Fs == -5):find(lags/Fs == -1),:));
    shadederrbar(lags/Fs, a, nanmean(corr_an{z,4}-corr_an{z,3},2), 'k'); hold on % SHUFFLE
    
    a = corr_an{z,1}; a = a - nanmean(a(find(lags/Fs == -5):find(lags/Fs == -1),:));
    shadederrbar(lags/Fs, nanmean(a,2), SEM(corr_an{z,1},2), clr{z}); % EARLY 
    
    a = corr_an_late{z,1}; a = a - nanmean(a(find(lags/Fs == -5):find(lags/Fs == -1),:));
    shadederrbar(lags/Fs, nanmean(a,2), SEM(corr_an_late{z,1},2), 'k'); % LATE
    
    xlabel('Lag (s)'); xlim([-0.5 0.5]);
    ylabel('Correlation Coefficient'); ylim([-1 0.5]); yticks([-1:0.25:1]);
    title(sprintf('%s (n = %d mice)',lbl{z}, nAn)); axis('square');
end
movegui(gcf,'center');

%% PLOT STATS VAL
a = 0.8; b = 1.2; r1 = a + (b-a).*rand(length(uni),1); % random values for jitter
a = 1.8; b = 2.2; r2 = a + (b-a).*rand(length(uni),1);
pl1 = min_val; pl2 = min_val_late; % plotting 1, plotting 2

fig = figure; fig.Position(3) = 1375; 
for z = 1:3
    subplot(1,3,z); hold on
    plot(r1, pl1(:,z), '.k', 'MarkerSize', 20);
    errorbar(1, nanmean(pl1(:,z)), SEM(pl1(:,z),1), '.r', 'MarkerSize', 20);
    plot(r2, pl2(:,z), '.k', 'MarkerSize', 20);
    errorbar(2, nanmean(pl2(:,z)), SEM(pl2(:,z),1), '.r', 'MarkerSize', 20);
    xlim([0.5 2.5]); xticks([1 2]); xticklabels({'Early','Late'}); 
    ylabel('Max Coefficient'); ylim([-1 0]); yticks([-1:0.25:0]);
    [~,p] = ttest(pl1(:,z),pl2(:,z));
    title(sprintf('%s (ttest: %1.2f)',lbl{z},p)); axis square
end
movegui(gcf,'center');

%% PLOT STATS LAG
a = 0.8; b = 1.2; r1 = a + (b-a).*rand(length(uni),1); % random values for jitter
a = 1.8; b = 2.2; r2 = a + (b-a).*rand(length(uni),1);
pl1 = min_lag.*1000; pl2 = min_lag_late.*1000; % plotting 1, plotting 2

fig = figure; fig.Position(3) = 1375; 
for z = 1:3
    subplot(1,3,z); hold on
    plot(r1, pl1(:,z), '.k', 'MarkerSize', 20);
    errorbar(1, nanmean(pl1(:,z)), SEM(pl1(:,z),1), '.r', 'MarkerSize', 20);
    plot(r2, pl2(:,z), '.k', 'MarkerSize', 20);
    errorbar(2, nanmean(pl2(:,z)), SEM(pl2(:,z),1), '.r', 'MarkerSize', 20);
    xlim([0.5 2.5]); xticks([1 2]); xticklabels({'Early','Late'}); 
    ylabel('Latency (ms)'); ylim([-300 0]); yticks([-500:100:0]);
    [~,p] = ttest(pl1(:,z),pl2(:,z));
    title(sprintf('%s (ttest: %1.3f)',lbl{z},p)); axis square
end
movegui(gcf,'center');