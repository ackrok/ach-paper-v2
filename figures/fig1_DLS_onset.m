a = [1:length(modACh)]; rmv = [3 17 20 22 23 30 36:39 41:43 48 50]; a(rmv) = [];
beh = modACh(a); % Extract recordings with reward
[align_full, time, ev] = plot_fp2event(beh,[-4 1],0); % Align photometry to events
align_adj = align_full;
for x = 1:length(beh); y = 1;
    align_adj{x,y} = align_adj{x,y} - nanmean(beh(x).FP{y});
    align_adj{x,y} = align_adj{x,y} - nanmean(align_adj{x,y}(find(time == -4):find(time == -2),:));
end

% n = X mice
tmp = {}; for x = 1:length(beh); tmp{x} = strtok(beh(x).rec,'-'); end
uni = unique(tmp); nAn = length(uni);
align_an = cell(length(uni),1); align_an_avg = [];
for x = 1:nAn
    ii = find(strcmp(tmp,uni{x}));
    align_an_avg(:,x) = nanmean([align_adj{ii,1}],2);
    align_an{x} = [align_adj{ii,1}];
end

%% Population AVERAGE
fig = figure; fig.Position(3) = 1000;
subplot(1,2,1); hold on
shadederrbar(time, nanmean(align_an_avg,2), SEM(align_an_avg,2), 'g');
ylabel('ACh3.0 fluorescence'); ylim([-2 8]); yticks([-2:2:8]);
yyaxis right; shadederrbar(vel_time, nanmean(vel_an_avg,2), SEM(vel_an_avg,2), 'k');
ylabel('Velocity (cm/s)'); ylim([-4 16]); yticks([-4:4:16]);
xlabel('Latency to Onset (s)'); xlim([-1 1]); xticks([-1:0.5:1]);
title(sprintf('ACh3.0 DLS Population Average (n = %d mice)',size(align_an_avg,2)));

% Population HEATMAP
tmp_heat = align_an_avg([find(time == -1):find(time == 1)],:);
m = max(tmp_heat); [a,b] = sort(m);
tmp_heat = tmp_heat(:,b); % Rank by maximum

subplot(1,2,2);
h = heatmap(tmp_heat');
h.Title = 'ACh3.0 DLS Population Average';
h.XLabel = 'Latency to Onset (s)'; h.YLabel = 'Animal';
h.GridVisible = 'off';
h.Colormap = parula; h.ColorLimits = [-2 14];

%% PLOT ALL ANIMALS to determine example
figure; 
for x = 1:length(align_an) 
    sp(x) = subplot(3,5,x); hold on
    plot(time, align_an{x}, 'Color', [0.05 0.75 0.45 0.5]);
    shadederrbar(time, nanmean(align_an{x},2), SEM(align_an{x},2), 'k');
    xlabel('Latency to Onset (s)'); xlim([-1 1]); xticks([-1:0.5:1]);
    ylabel('ACh3.0 fluorescence'); %ylim([-7 17]); yticks([-5:5:15]);
    title(sprintf('%s (n = %d bouts)',uni{x},size(align_an{x},2)));
end

%% Example session INDIVIDUAL TRIALS
x = 3; % example

idx = [1:size(align_an{x},2)]';
idx(any(abs(vel_an{x}(find(vel_time == -1):find(vel_time == -0.5),:)) > 1)) = [];

fig = figure; fig.Position(3) = 1000; movegui(gcf,'center');
subplot(1,3,1); hold on
plot(time, align_an{x}(:,idx), 'Color', [0.05 0.75 0.45 0.5]);
shadederrbar(time, nanmean(align_an{x},2), SEM(align_an{x},2), 'k');
xlabel('Latency to Onset (s)'); xlim([-1 1]); xticks([-1:0.5:1]);
ylabel('ACh3.0 fluorescence'); %ylim([-7 17]); yticks([-5:5:15]);
title(sprintf('%s (n = %d bouts)',uni{x},length(idx)));
axis('square');

subplot(1,3,2); hold on
plot(vel_time, vel_an{x}(:,idx), 'Color', [0 0 0 0.1]);
shadederrbar(vel_time, nanmean(vel_an{x},2), SEM(vel_an{x},2), 'k');
xlabel('Latency to Onset (s)'); xlim([-1 1]); xticks([-1:0.5:1]);
ylabel('Velocity (cm/s)');
axis('square');

subplot(1,3,3);
tmp_heat = align_an{x}(find(time == -1):find(time == 1),idx);
m = max(align_an{x}(find(time == 0):find(time == 0.4),idx)); 
[~,b] = sort(m);
tmp_heat = tmp_heat(:,b); % Rank by maximum

h = heatmap(tmp_heat');
h.Title = 'ACh3.0 DLS Population Average';
h.XLabel = 'Latency to Onset (s)'; h.YLabel = 'Animal';
h.GridVisible = 'off';
h.Colormap = parula; h.ColorLimits = [-5 20];