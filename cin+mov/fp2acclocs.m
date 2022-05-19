a = [1:length(modACh)]; rmv = [3 17 20 22 23 30 36:39 41:43 48 50]; a(rmv) = [];
beh = modACh(a); % Extract recordings with reward
[align_full, time, ev] = plot_fp2event(beh,[-4 1],0); % Align photometry to events
%
align_adj = align_full;
for x = 1:length(beh); y = 1;
    align_adj{x,y} = align_adj{x,y} - nanmean(beh(x).FP{y});
    align_adj{x,y} = align_adj{x,y} - nanmean(align_adj{x,y}(find(time == -4):find(time == -1),:));
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
xlabel('Latency to Acc Peak (s)'); xlim([-1 1]); xticks([-1:0.5:1]);
title(sprintf('ACh3.0 (n = %d mice)',size(align_an_avg,2)));

% Population HEATMAP
tmp_heat = align_an_avg([find(time == -1):find(time == 1)],:);
m = max(tmp_heat); [~,m_sort] = sort(m);
tmp_heat = tmp_heat(:,m_sort); % Rank by maximum
subplot(1,2,2);
h = heatmap(tmp_heat');
h.Title = 'Average';
h.XLabel = 'Latency to Acc Peak (s)'; h.YLabel = 'Animal';
h.GridVisible = 'off';
h.Colormap = jet; h.ColorLimits = [-6 12];

%% SHUFFLE
align_shuff_full = cell(length(beh),3);
nShuff = 50;
for x = 1:length(beh)
    acc = getAcc(beh(x).vel); % Extract acceleration signal
    [~,locs] = findpeaks(acc,'MinPeakProminence',0.5,'MinPeakDistance',0.5); % Location of peaks, using findpeaks function
    ev = beh(x).time(locs); % Convert peak locations to seconds
    clc
    sig = beh(x).FP{1}; Fs = 50; % Signal that will be aligned to event times
    for z = 1:nShuff
        ev_rand = randperm(length(sig),length(ev)); % random sampling without replacement
        mat = getSTA(sig, ev_rand, Fs, [-4 1]);
        prc = prctile(mat,[5 50 95],2); %5th, 50th, 95th percentile of STA
        align_shuff_full{x,1}(:,z) = prc(:,1);
        align_shuff_full{x,2}(:,z) = prc(:,2);
        align_shuff_full{x,3}(:,z) = prc(:,3);
    end
end

%% n = X mice
tmp = {}; for x = 1:length(beh); tmp{x} = strtok(beh(x).rec,'-'); end
uni = unique(tmp); nAn = length(uni);
align_shuff = cell(length(uni),3); 
align_lowupp = cell(length(uni),2); % first column is 5th prc and second column is 95th prc
for x = 1:nAn
    ii = find(strcmp(tmp,uni{x}));
    for z = 1:3
        align_shuff{x,z} = [align_shuff_full{ii,z}];
    end
    align_lowupp{x,1} = nanmean(align_shuff{x,1} - align_shuff{x,2},2);
    align_lowupp{x,2} = nanmean(align_shuff{x,3} - align_shuff{x,2},2);
end

%% PLOT proportion > 95% CI
above95 = []; below5 = [];
for x = 1:nAn
    a = align_an_avg(:,x); b = align_lowupp{x,2}; c = align_lowupp{x,1};
    above95 = [above95, a > b]; %binary vector where CCG passed 95% confidence interval
    below5 = [below5, a < c]; %binary vector where CCG below 5% confidence interval
end
figure; hold on % figure; hold on
ds = 2;
a = 100*sum(above95,2)/size(above95,2); 
[prop_m,prop_t] = max(a); prop_t = time(prop_t); 
a = a(1:ds:end);
b = -100*sum(below5,2)/size(below5,2); b = b(1:ds:end);
bar(time(1:ds:end), a,'FaceColor','b','FaceAlpha',0.5);
bar(time(1:ds:end), b,'FaceColor','b','FaceAlpha',0.5);
xlabel('Lag (s)'); ylabel('Prop of Pair CCG > 95% CI (%)'); 
ylim([-50 100]); xlim([-1 1]);
title(sprintf('max %1.1f @ %d ms (n = %d units)',prop_m, round(1000*prop_t), size(above95,2)))
axis('square');
