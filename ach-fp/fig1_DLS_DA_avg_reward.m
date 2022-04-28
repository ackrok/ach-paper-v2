good_rew = [10:16,20,22:26,28:29,32:35,40:42,44:46];
beh = modAChDA(good_rew);

%%
y = 2;
[align_full, time, ev] = plot_fp2event(beh,[-6 2],0); % Align photometry to events
align_adj = align_full;
for x = 1:length(beh)
    align_adj{x,y} = align_adj{x,y} - nanmean(beh(x).FP{y}); 
    align_adj{x,y} = align_adj{x,y} - nanmean(align_adj{x,y}(find(time == -6):find(time == -1),:),1);
end

%n = X mice
tmp = {}; for x = 1:length(beh); tmp{x} = strtok(beh(x).rec,'-'); end
uni = unique(tmp); nAn = length(uni);
align_an = cell(length(uni),1); align_an_avg = [];
for x = 1:nAn
    ii = find(strcmp(tmp,uni{x}));
    align_an_avg(:,x) = nanmean([align_adj{ii,y}],2);
    align_an{x} = [align_adj{ii,y}];
end

%% Population AVERAGE + HEATMAP
fig = figure; fig.Position(3) = 1000;
subplot(1,2,1); hold on
shadederrbar(time, nanmean(align_an_avg,2), SEM(align_an_avg,2), 'm');
% plot(time, align_an_avg, 'Color', [0 0 0 0.1]);
xlabel('Latency to reward (s)'); xlim([-1 2]); xticks([-1:0.5:2]);
ylabel('DA fluorescence'); ylim([-2 10]); yticks([-2:2:10]);
title(sprintf('DA DLS Population Average (n = %d mice)',size(align_an_avg,2)));
% axis('square');

%Population HEATMAP
m = max(align_an_avg);
[a,b] = sort(m);
subplot(1,2,2);
h = heatmap(align_an_avg([find(time == -1):find(time == 2)],flip(b))');
h.Title = 'DA DLS Population Average';
h.XLabel = 'Latency to reward (s)'; h.YLabel = 'Animal';
h.GridVisible = 'off';
h.Colormap = parula; h.ColorLimits = [-2 12];

%Population HEATMAP
% m = max(align_an_avg);
% [a,b] = sort(m);
% subplot(1,3,3);
% h = heatmap(align_an_avg([find(time == -1):find(time == 2)],b)');
% h.Title = 'ACh3.0 DLS Population Average';
% h.XLabel = 'Latency to reward (s)'; h.YLabel = 'Animal';
% h.GridVisible = 'off';
% h.Colormap = jet; h.ColorLimits = [-6 10];

%% Example session INDIVIDUAL TRIALS
% for x = 1:length(align_an)
    
x = 1; % example

idx = [1:size(align_an{x},2)]';
a = align_an{x}(1,:); % first row
idx = idx(~isnan(a)); % remove nans
if length(idx) > 100; idx = datasample(idx,100,'Replace',false); end % choose 100 random trials
m = min(align_an{x}(find(time == 0):find(time == 1),idx)); % minima of 100 trials
[~,b] = sort(m); % sort by reward minima

fig = figure; %fig.Position(3) = 1000;
% subplot(4,3,x);
hold on
plot(time, align_an{x}(:,idx), 'Color', [0.05 0.75 0.45 0.5]);
shadederrbar(time, nanmean(align_an{x},2), SEM(align_an{x},2), 'k');
xlabel('Latency to reward (s)'); xlim([-1 2]); xticks([-1:0.5:2]);
ylabel('DA fluorescence'); % ylim([-7 17]); yticks([-5:5:15]);
title('DA DLS Example (n = 100 trials)');
% title(sprintf('%s',uni{x}));
% end

% figure;
% tmp = align_an{x}([find(time == -1):find(time == 2)],flip(idx(b)))';
% h2 = heatmap(tmp);
% h2.Title = 'ACh3.0 DLS Example';
% h2.XLabel = 'Latency to reward (s)'; h2.YLabel = 'Trial';
% h2.GridVisible = 'off';
% h2.Colormap = jet; h2.ColorLimits = [-5 30];

%% Licking
Fs = 50; bin = 0.1;
align_lick = cell(length(beh),1);
for x = 1:length(beh)
    lick = beh(x).lick(:)/beh(x).Fs; % Licks, in seconds
    lickms = lick*1000;
    a2 = diff(lickms);
    a3 = find(a2 < 2);
    lick(a3) = nan;
    peth = getClusterPETH(lick, ev{x}(~isnan(ev{x})), bin, [-1 2]);
    peth_base = getClusterPETH(lick, ev{x}(~isnan(ev{x})), bin, [-3 -2]);
    align_lick{x} = peth.cts{1}./(bin);
    % align_lick{x} = align_lick{x} - nanmean(peth_base.cts{1}); 
end
tmp = {}; for x = 1:length(beh); tmp{x} = strtok(beh(x).rec,'-'); end
uni = unique(tmp); nAn = length(uni);
lick_an = cell(length(uni),1); lick_an_avg = [];
for x = 1:nAn
    ii = find(strcmp(tmp,uni{x}));
    lick_an{x} = [align_lick{ii,1}];
    lick_an_avg(:,x) = nanmean([align_lick{ii,1}],2);
end

%% LICK population AVERAGE
sm = 1; 
fig = figure; fig.Position(3) = 1000;
subplot(1,2,1);
a = nanmean(lick_an_avg,2);
% a(1:37) = a(1:37) - nanmean(a(1:37));
% a = movmean(a,sm);
% a(38:50) = 0;
shadederrbar(peth.time, a, movmean(SEM(lick_an_avg,2),sm), 'k');
xlabel('Latency to reward (s)'); xlim([-1 2]); xticks([-1:0.5:2]);
ylabel('Lick (Hz)'); % ylim([0 1]); yticks([0:0.2:1]);
title(sprintf('Lick Population Average (n = %d mice)',size(lick_an_avg,2)));

% LICK population HEATMAP
subplot(1,2,2);
m = min(align_an_avg);
[a,b] = sort(m);
h = heatmap(lick_an_avg(:,b)');
h.Title = 'Lick Population Average';
h.XLabel = 'Latency to reward (s)'; h.YLabel = 'Animal';
h.GridVisible = 'off';
h.Colormap = parula; h.ColorLimits = [0 20];

%% LICK example session RASTER
% use x from example above
x = 3;
zz = sort(randi(size(lick_an{x},2),[100,1]));
figure; hold on
plot([0 0],[0 100],'-k');
for ii = 1:length(zz)
    plot(peth.time, ii.*(lick_an{x}(:,ii) > 0), '.k');
end
xlabel('Latency to reward (s)'); xlim([-1 2]); xticks([-1:0.5:2]);
ylabel('Trial Number');
title('ACh3.0 DLS Example (n = 100 trials)');

%% TRACE EXAMPLE
x = 15;
figure; hold on
plot(beh(x).time, beh(x).FP{1} - nanmean(beh(x).FP{1}), 'g');
plot(beh(x).time, getAcc(beh(x).vel), 'k');
stem(beh(x).reward/beh(x).Fs, ones(length(beh(x).reward),1), 'b');
plot(beh(x).lick/beh(x).Fs, 2*ones(length(beh(x).lick),1), '.k');
xlim([340 355]);
xlim([522.6 544.6]);
xlabel('Time (s)'); ylabel('ACh3.0 Fluorescence');
title(sprintf('%s',beh(x).rec));
legend({'ACh3.0','acc','rew','lick'})

%% How many are pause vs burst-pause?
figure;
for x = 1:length(align_an)
    subplot(5,3,x); hold on
	
    tmp = align_an{x};
    time_seg = time([find(time == -1):find(time == 2)])';
    align_seg = tmp([find(time == -1):find(time == 2)],:);
    base = tmp([find(time == -6):find(time == -2)],:);
    sig = nanstd(nanmean(base,2)); % 1 standard deviation

    subplot(5,3,x); hold on
    shadederrbar(time_seg, nanmean(align_seg,2), SEM(align_seg,2), 'g');
    shadederrbar(time_seg, zeros(length(time_seg),1), 2*sig.*ones(length(time_seg),1), 'k');
    title(sprintf('%s',uni{x}));
end