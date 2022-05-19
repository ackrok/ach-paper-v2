%% Extract recordings
good_rew = [10:16,20,22:26,28:29,32:35,40:42,44:46];
beh = modAChDA(good_rew);

% a = [10:16,20,22:30,36:37,39:43,44:49];
% beh = modACh(a); % Extract recordings with reward

%% Align data to reward
[align_full, time, ev] = plot_fp2event(beh,[-6 2],0); % Align photometry to events
Fs = 50;

align_adj = align_full;
for x = 1:length(beh)
    for y = 1:2
    align_adj{x,y} = align_adj{x,y} - nanmean(beh(x).FP{y}); 
    align_adj{x,y} = align_adj{x,y} - nanmean(align_adj{x,y}(find(time == -6):find(time == -1),:),1);
    end
end

% n = X mice
tmp = {}; for x = 1:length(beh); tmp{x} = strtok(beh(x).rec,'-'); end
uni = unique(tmp); nAn = length(uni);
align_an = cell(length(uni),2); 
align_avg = cell(1,2); align_norm = cell(1,2);
lag = []; val = []; 
fwhm = []; win = 0.3*Fs; % set window for pause duration analysis, usually 0.3*Fs (resulting window is 2*win)
for x = 1:nAn
    ii = find(strcmp(tmp,uni{x}));
    for y = 1:2
        align_avg{y}(:,x) = nanmean([align_adj{ii,y}],2);
        align_an{x,y} = [align_adj{ii,y}];

        vec = nanmean(align_an{x,y},2);
        bb = nanmean(vec([find(time==-6):find(time==-1)],:),1);
        
        if y == 1 
            mm = min(vec([find(time==0):find(time==1)],:));
            norm = (vec - mm)./(bb - mm);
        elseif y == 2
            mm = max(vec([find(time==0):find(time==1)],:));
            norm = (vec - bb)./(mm - bb); end
        align_norm{y}(:,x) = norm; % normalize(align_an_avg{1,y}(:,x),'range'); % Normalize so range is [0 1]
    end

    y = 1;
    r = [0.1 0.3]; Fs = 50; r = (r/(1/Fs) + find(time==0)); % range for ACh peak
    [a,b] = max(align_an{x,y}(r(1):r(2),:)); % find local maximum within range
    c = time(b + r(1) - 1); % convert index to seconds
    c(b == 1) = nan; c(isnan(a)) = nan; % remove values that are not local maxima
    lag(x,1) = nanmean(c);
    val(x,1) = nanmean(a);
    
    y = 2;
    r = [0.1 0.5]; r = (r/(1/Fs) + find(time==0)); % range for DA peak
    [a,b] = max(align_an{x,y}(r(1):r(2),:)); % find local maximum within range
    c = time(b + r(1) - 1); c(isnan(a)) = nan; % convert index to seconds
    lag(x,2) = nanmean(c);
    val(x,2) = nanmean(a);
    
    y = 1;
    r = [0.2 0.7]; r = (r/(1/Fs) + find(time==0)); % range for ACh trough
    [a,b] = min(align_an{x,y}(r(1):r(2),:)); % find local minima within range
    c = time(b + r(1) - 1); c(isnan(a)) = nan; % convert index to seconds
    lag(x,3) = nanmean(c);
    val(x,3) = nanmean(a);
    
    for y = 1 % duration of ACh trough
        tmp_dur_pause = nan(size(align_an{x,y},2),1);
        base = align_an{x,y}(find(time==-6):find(time==-1),:);
        thres = 1.5*(nanstd(base(:)));
        for z = 1:size(align_an{x,y},2)
            signal = align_an{x,y}(:,z);
            [mm,ii] = min(signal(find(time==0):find(time==1)));
            if abs(mm) < thres; continue; end % continue if pause does not pass threshold
            ii = ii + find(time==0) -1;
            halfMax = 0.5*mm; % amplitude at half-max
            a = [signal(ii-win : ii)]; a = flipud(a); % segment preceding idx of maximum deflection
            b = find(a > 0, 1, 'first') - 1;
            if ~isempty(b); a = a(1:b); end % retain only values above 0
            c = []; [~,c(1)] = min(abs(a-halfMax)); c(1) = c(1)-1; % idx at half max
            a = [signal(ii : ii+win)]; % segment following idx of maximum deflection
            b = find(a > 0, 1, 'first') - 1;
            if ~isempty(b); a = a(1:b); end % retain only values above 0
            [~,c(2)] = min(abs(a-halfMax)); % idx at half max
            tmp_dur_pause(z) = sum(c);
        end
        fwhm(x,y) = nanmean(tmp_dur_pause);
    end
end
lag([12 13],1) = nan;

%% STATISTICS
y = 2;
[a,b] = max(align_avg{y});
c = SEM(align_avg{y},2);
d = lag(:,2).*1000;
fprintf('DA DLS peak amplitude: %1.2f +/- %1.2f dF/F , %1.2f +/- %1.2f ms from reward (n=%d mice, %d recs) \n',nanmean(a),c(round(nanmean(b))),nanmean(d),SEM(d,1),size(align_avg{y},2),length(beh));

y = 1;
[a,b] = min(align_avg{y});
c = SEM(align_avg{y},2);
d = lag(:,3).*1000;
fprintf('ACh DLS pause amplitude: %1.2f +/- %1.2f dF/F , %1.2f +/- %1.2f ms from reward (n = %d mice, %d recs) \n',nanmean(a),c(round(nanmean(b))),nanmean(d),SEM(d,1),size(align_avg{y},2),length(beh));
fprintf('ACh pause duration: %1.2f +/- %1.2f ms \n',nanmean(fwhm*(1000/50)),SEM(fwhm*(1000/50),1));

[a,b] = max(align_avg{y}); 
a = a(b > find(time == 0)); b = b(b > find(time == 0)); % restrict to animals that have burst-pause response
c = SEM(align_avg{y},2);
d = lag([1:length(b)],1).*1000;
fprintf('ACh DLS peak amplitude: %1.2f +/- %1.2f dF/F , %1.2f +/- %1.2f ms from reward (n = %d mice) \n',nanmean(a),c(round(nanmean(b))),nanmean(d),SEM(d,1),length(b));

a = (lag(:,2)-lag(:,1)).*1000;
fprintf('ACh peak to DA peak: %1.2f +/- %1.2f ms \n',nanmean(a),SEM(a,1));

a = (lag(:,3)-lag(:,2)).*1000;
fprintf('DA peak to ACh pause: %1.2f +/- %1.2f ms \n',nanmean(a),SEM(a,1));

[p,~,stats] = kruskalwallis(lag,[],'off');
c = multcompare(stats,'display','off');
fprintf('latency to reward kruskalwallis p = %1.3f, AChpeak/DApeak p = %1.3f, DApeak/AChtrough p = %1.3f \n',p,c(1,6),c(3,6));

%% Population AVERAGE + HEATMAP
y = 2;
fig = figure; fig.Position([3 4]) = [850 800];
subplot(2,2,1); hold on
shadederrbar(time, nanmean(align_avg{y},2), SEM(align_avg{y},2), 'm');
% plot(time, align_avg{y}, 'Color', [0 0 0 0.1]);
xlabel('Latency to reward (s)'); xlim([-1 2]); xticks([-1:0.5:2]);
ylabel('DA fluorescence'); ylim([-2 12]); yticks([-2:2:12]);
[a,b] = max(align_avg{y});
c = SEM(align_avg{y},2);
title(sprintf('DA DLS: %1.2f +/- %1.2f max (n = %d mice)',nanmean(a),c(round(nanmean(b))),size(align_avg{y},2)));
axis('square');

%Population HEATMAP
m = max(align_avg{y});
[sortM,sortIdx] = sort(m);
subplot(2,2,2);
h = heatmap(align_avg{y}([find(time == -1):find(time == 2)],flip(sortIdx))');
h.Title = 'DA DLS Population Average';
h.XLabel = 'Latency to reward (s)'; h.YLabel = 'Animal';
h.GridVisible = 'off';
h.Colormap = parula; h.ColorLimits = [-2 14];

y = 1; % ACh
subplot(2,2,3); hold on
shadederrbar(time, nanmean(align_avg{y},2), SEM(align_avg{y},2), 'g');
% plot(time, align_avg{y}, 'Color', [0 0 0 0.1]);
xlabel('Latency to reward (s)'); xlim([-1 2]); xticks([-1:0.5:2]);
ylabel('ACh fluorescence'); ylim([-6 4]); yticks([-6:2:4]);
[a,b] = min(align_avg{y});
c = SEM(align_avg{y},2);
title(sprintf('ACh DLS: %1.2f +/- %1.2f min (n = %d mice)',nanmean(a),c(round(nanmean(b))),size(align_avg{y},2)));;
axis('square');

%Population HEATMAP
subplot(2,2,4);
h2 = heatmap(align_avg{y}([find(time == -1):find(time == 2)],flip(sortIdx))');
h2.Title = 'ACh DLS Population Average';
h2.XLabel = 'Latency to reward (s)'; h2.YLabel = 'Animal';
h2.GridVisible = 'off';
h2.Colormap = parula; h2.ColorLimits = [-6 10];

movegui(gcf,'center')

%% OVERLAY POPULATION PLOT
%PLOT: %dF/F 
fig = figure; fig.Position(3) = 1375;
subplot(1,3,1);
plot([0 0],[-5 11],'--k'); hold on
shadederrbar(t, nanmean(align_avg{2},2), SEM(align_avg{2},2), 'm'); hold on
shadederrbar(t, nanmean(align_avg{1},2), SEM(align_avg{1},2), 'g');
ylabel('FP (% dF/F)'); xlabel('Latency to reward delivery (s)');
xlim([-1 2]);
title(sprintf('reward: ACh3.0 + rDA1m (n = %d mice)',size(align_an,1))); axis('square');
 
%PLOT: %dF/F normalized range [0 1]
subplot(1,3,2);
plot([0 0],[-1.1 1.1],'--k'); hold on
shadederrbar(t, nanmean(align_norm{2},2), SEM(align_norm{2},2), 'm'); hold on
shadederrbar(t, nanmean(align_norm{1}-1,2), SEM(align_norm{1},2), 'g');
ylabel('FP (a.u.)'); xlabel('Latency to reward delivery (s)');
xlim([-0.2 0.8]); xticks([0 0.5]);
ylim([-1.1 1.1]); yticks([-1:1]);
title(sprintf('reward: ACh3.0 + rDA1m (n = %d mice)',size(align_an,1)));  axis('square');

subplot(1,3,3); hold on;
a = lag'.*1000;
plot(a,'--','Color',[0 0 0 0.2]);
errorbar(nanmean(a,2),SEM(a,2),'.-k','MarkerSize',20);
xlim([0.5 3.5]); xticks([1:3]); xticklabels({'ACh peak','DA peak','ACh trough'});
ylabel('latency to reward (ms)'); ylim([0 700]);
% p = []; p(1) = signrank(a(:,1),a(:,2)); p(2) = signrank(a(:,2),a(:,3));
p = kruskalwallis(a',[],'off');
title(sprintf('kruskalwallis: p = %1.3f',p)); axis('square');
movegui(gcf,'center')

%% Example session INDIVIDUAL TRIALS
x = 1; % example AK189
% x = 8; % example JM003R
nTraces = 20; % how many traces to plot
fig = figure; fig.Position(3) = 1000; clr = {'g','m'};

idx = [1:size(align_an{x,1},2)]';
a = align_an{x,1}(1,:); % first row
idx = idx(~isnan(a)); % remove nans
if length(idx) > nTraces; idx = datasample(idx,nTraces,'Replace',false); end % choose N random trials

for y = 1:2
    % m = min(align_an{x,y}(find(time == 0):find(time == 1),idx)); % minima of N trials
    % [~,b] = sort(m); % sort by reward minima

    subplot(1,2,y); hold on
    plot(time, align_an{x,y}(:,idx), 'Color', [0 0 0 0.1]);
    shadederrbar(time, nanmean(align_an{x,y},2), SEM(align_an{x,y},2), clr{y});
    xlabel('Latency to reward (s)'); xlim([-1 2]); xticks([-1:0.5:2]);
    ylabel('fluorescence'); % ylim([-7 17]); yticks([-5:5:15]);
    title(sprintf('example %s (n = %d trials)',uni{x},nTraces));
    % title(sprintf('%s',uni{x}));
end
movegui(gcf,'center')
% figure;
% tmp = align_an{x}([find(time == -1):find(time == 2)],flip(idx(b)))';
% h2 = heatmap(tmp);
% h2.Title = 'ACh3.0 DLS Example';
% h2.XLabel = 'Latency to reward (s)'; h2.YLabel = 'Trial';
% h2.GridVisible = 'off';
% h2.Colormap = jet; h2.ColorLimits = [-5 30];