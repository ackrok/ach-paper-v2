%% Extract recordings
load('C:\Users\Anya\Desktop\FP_LOCAL\AK189-197_cannula+rew_v4.mat')
beh = s(1).s; % CHANGE

%% Licking process
[align_full, time, ev] = plot_fp2event(beh,[-6 2],0); % Align to events
Fs = 50; 
align_lick = cell(length(beh),1); lick_first = align_lick;
for x = 1:length(beh)
    lick = beh(x).lick(:)/beh(x).Fs; % Licks, in seconds
    lick_repeat = [diff(lick.*1000) > 50]; % Identify licks that are <50ms after previous lick
    lick_sub = lick; lick_sub(1) = [];
    lick = [lick(1); lick_sub(lick_repeat)];
    peth = getClusterPETH(lick, ev{x}(~isnan(ev{x})), 1/Fs, [-1 2]);
    pt = peth.time;
    % peth_base = getClusterPETH(lick, ev{x}(~isnan(ev{x})), 1/Fs, [-3 -1]);
    align_lick{x} = peth.cts{1}./(1/Fs);
    % align_lick{x} = align_lick{x} - nanmean(peth_base.cts{1}); 
    % align_lick{x} = align_lick{x} - nanmean(nanmean(align_lick{x}(1:50,:),2));
    
    bin = 1/1000; window = [0 0.5];
    peth = getClusterPETH(lick, ev{x}(~isnan(ev{x})), bin, window); % PETH: lick aligned to reward in 1 ms bins
    cts = peth.cts{1}; % Lick counts in 1ms bins for each reward trial
    [~, lick_first{x}] = max(cts~=0, [], 1); % Find first non-zero index for each trial
end
tmp = {}; for x = 1:length(beh); tmp{x} = strtok(beh(x).rec,'-'); end
uni = unique(tmp); nAn = length(uni);
lick_an = cell(length(uni),1); lick_an_avg = []; lick_0_an = lick_an;
for x = 1:nAn
    ii = find(strcmp(tmp,uni{x}));
    lick_an{x} = [align_lick{ii,1}];
    lick_an_avg(:,x) = nanmean([align_lick{ii,1}],2);
    lick_0_an{x} = [lick_first{ii,1}];
end

%% LICK population AVERAGE
fig = figure; fig.Position(3) = 1000;
subplot(1,2,1); sm = 10;
shadederrbar(pt, movmean(nanmean(lick_an_avg,2),sm), movmean(SEM(lick_an_avg,2),sm), 'k');
xlabel('Latency to reward (s)'); xlim([-1 2]); xticks([-1:0.5:2]);
ylabel('Lick (Hz)'); % ylim([0 1]); yticks([0:0.2:1]);
title(sprintf('Lick Population Average (n = %d mice)',size(lick_an_avg,2)));

% LICK population HEATMAP
subplot(1,2,2);
% m = min(lick_an_avg);
% [a,b] = sort(m);
% b = randperm(size(lick_an_avg,2));
h = heatmap(movmean(lick_an_avg',10,2));
h.Title = 'Lick Population Average';
h.XLabel = 'Latency to reward (s)'; h.YLabel = 'Animal';
h.GridVisible = 'off';
h.Colormap = gray; % caxis([0 1.5])