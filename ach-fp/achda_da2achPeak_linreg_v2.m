%% run achda_da2achpeak

amp = cell(size(ach2achpeak,1),1); % initialize matrix for amplitude of ach peak
lag = cell(size(ach2achpeak,1),1); % latency of ach peak to check that it is ~0ms

for x = 1:size(ach2achpeak,1)
    [m, ii] = max(ach2achpeak{x,1});
    amp{x,1} = m(:);
    lag{x,1} = sta_time(ii(:));
    
    m1 = min(da2achpeak{x,1}(find(sta_time == -0.26):find(sta_time == 0),:)); % find MIN within range preceding ACh peak
    [m2, ii] = max(da2achpeak{x,1}(find(sta_time == 0):find(sta_time == 0.26),:)); % find MAX within range following ACh peak
    amp{x,2} = m2(:) - m1(:); % difference in amplitude
    lag{x,2} = sta_time(ii(:) + find(sta_time == 0) - 1);
end

% N = X mice: DA to ACh pause
tmp = {}; for x = 1:length(beh); tmp{x} = strtok(beh(x).rec,'-'); end
uni = unique(tmp); nAn = length(uni);
amp_an = cell(nAn,2); lag_an = cell(nAn,2);
for x = 1:nAn; for y = 1:2
    ii = find(strcmp(tmp,uni{x}));
    for z = 1:length(ii)
        amp_an{x,y} = [amp_an{x,y}; amp{ii(z),y}];
        lag_an{x,y} = [lag_an{x,y}; lag{ii(z),y}];
    end
    end
end

% Slope, p-value
stats = {'r','p','slope'}; stats_amp = []; stats_lag = [];
fit_amp = []; fit_lag = [];
for x = 1:nAn
    mdl_amp = fitlm(amp_an{x,1}, amp_an{x,2}); % Linear regression: DA ~ ACh
    mdl_lag = fitlm(amp_an{x,1}, lag_an{x,2}); 

    y = 1; stats_amp(x,y) = mdl_amp.Rsquared.Adjusted; stats_lag(x,y) = mdl_lag.Rsquared.Adjusted;
    y = 2; stats_amp(x,y) = mdl_amp.Coefficients{2,4}; stats_lag(x,y) = mdl_lag.Coefficients{2,4};
    y = 3; stats_amp(x,y) = mdl_amp.Coefficients{2,1}; stats_lag(x,y) = mdl_lag.Coefficients{2,1};
    fit_amp(:,x) = mdl_amp.Coefficients{1,1} + mdl_amp.Coefficients{2,1}.*[0:20];
    fit_lag(:,x) = mdl_lag.Coefficients{1,1} + mdl_lag.Coefficients{2,1}.*[0:20];
end

%% PLOT ALL
figure;
for x = 1:length(uni) 
    sp(x) = subplot(3,2,x); hold on
    plot(amp_an{x,1},amp_an{x,2},'.k');
    plot([0:20], fit_amp(:,x), 'r', 'DisplayName', 'lin fit');
    xlabel('ACh peak amplitude'); ylabel('DA peak amplitude');
    title(sprintf('%s',uni{x})); axis('equal');
end
linkaxes(sp,'y'); linkaxes(sp,'x');

%% PLOT EXAMPLE
fig = figure; fig.Position(3) = 1375;

% x = 9; % JM004 control
% x = 4; % AK188 b2flox
subplot(1,3,1); hold on
x = 9; % JM004 control
z = randsample(length(amp_wt{x,1}),250);
plot(amp_wt{x,1}(z),amp_wt{x,2}(z),'.k', 'DisplayName', 'data');
plot([0:20], fit_amp(:,x), 'k', 'DisplayName', 'lin fit');
x = 4; % AK188 b2flox
z = randsample(length(amp_b2{x,1}),250);
plot(amp_b2{x,1}(z),amp_b2{x,2}(z),'.r', 'DisplayName', 'data');
plot([0:20], fit_amp(:,x), 'r', 'DisplayName', 'lin fit');
xlabel('ACh peak amplitude'); ylabel('DA peak amplitude');
title(sprintf('control: r = 0.011 | p = 0.055 | slope = 0.11 \n b2flox: r = 0.050 | p = 0.000 | slope = 0.21'))
% title(sprintf('%s: r = %1.3f | p = %1.3f | slope = %1.2f',uni{x},...
%     stats_amp(x,1),stats_amp(x,2),stats_amp(x,3))); 
axis('equal');
movegui(gcf,'center');

subplot(1,3,2); hold on
violinplot(stats_wt); ylim([-0.1 1.1])
xticklabels(stats);
title(sprintf('AVG: r = %1.3f | p = %1.3f | slope = %1.2f',...
    nanmean(stats_wt(:,1)),nanmean(stats_wt(:,2)),nanmean(stats_wt(:,3))));

subplot(1,3,3); hold on
violinplot(stats_b2); ylim([-0.1 1.1])
xticklabels(stats);
title(sprintf('b2flox: r = %1.3f | p = %1.3f | slope = %1.2f',...
    nanmean(stats_b2(:,1)),nanmean(stats_b2(:,2)),nanmean(stats_b2(:,3))));

%%
% stats_b2 = stats_amp; 
% amp_b2 = amp_an;
% fit_b2 = fit_amp;

stats_wt = stats_amp;
amp_wt = amp_an; lag_wt = lag_an;
fit_wt = fit_amp;

%%
p = [];
for y = 1:3; p(y) = ranksum(stats_wt(:,y), stats_b2(:,y)); end

%% PLOT ALL
figure;
for x = 1:length(lag_wt) 
    sp(x) = subplot(3,4,x); hold on
    violinplot(lag_wt{x,2});
    title(sprintf('%s',uni{x}));
end
linkaxes(sp,'y');
