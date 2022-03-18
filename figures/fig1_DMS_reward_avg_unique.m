a = [7:12,16];
beh = allDMS(a); % Extract recordings with reward
[align_full, t] = plot_fp2event(beh,[-6 2]); % Align photometry to events

%%
align_adj = align_full; align_an_avg = [];
for x = 1:length(beh)
    align_adj{x,1} = align_adj{x,1} - nanmean(beh(x).FP{1}); 
    align_adj{x,1} = align_adj{x,1} - nanmean(align_adj{x,1}(find(time == -6):find(time == -2),:),1);
    align_an_avg(:,x) = nanmean(align_adj{x,1},2);
end

%% Population AVERAGE + HEATMAP
fig = figure; fig.Position(3) = 1000;
subplot(2,2,1); hold on
shadederrbar(time, nanmean(align_an_avg,2), SEM(align_an_avg,2), 'g');
% plot(time, align_an_avg, 'Color', [0 0 0 0.1]);
xlabel('Latency to reward (s)'); xlim([-1 2]); xticks([-1:0.5:2]);
ylabel('ACh3.0 fluorescence'); ylim([-2 1]); yticks([-4:4]);
title(sprintf('ACh3.0 DMS Pop Average (n = %d recs)',size(align_an_avg,2)));
axis('square');

%Population HEATMAP
m = min(align_an_avg);
[a,b] = sort(m);
subplot(2,2,2);
h = heatmap(align_an_avg([find(time == -1):find(time == 2)],flip(b))');
h.Title = 'ACh3.0 DMS Pop Average';
h.XLabel = 'Latency to reward (s)'; h.YLabel = 'Animal';
h.GridVisible = 'off';
h.Colormap = parula; h.ColorLimits = [-3 1];

% EXAMPLE
x = 6; % example

idx = [1:size(align_adj{x},2)]';
a = align_adj{x}(1,:); % first row
idx = idx(~isnan(a)); % remove nans
% idx = datasample(idx,100,'Replace',false); % choose 100 random trials
m = min(align_adj{x}(find(time == 0):find(time == 1),idx)); % minima of 100 trials
[~,b] = sort(m); % sort by reward minima

% fig = figure; fig.Position(2) = 1000;
subplot(2,2,3); hold on
plot(time, align_adj{x}(:,idx), 'Color', [0.05 0.75 0.45 0.5]);
shadederrbar(time, nanmean(align_adj{x},2), SEM(align_adj{x},2), 'k');
xlabel('Latency to reward (s)'); xlim([-1 2]); xticks([-1:0.5:2]);
ylabel('ACh3.0 fluorescence'); % ylim([-7 17]); yticks([-5:5:15]);
title(sprintf('ACh3.0 DMS Example (n = %d trials)',length(idx)));
axis('square');

subplot(2,2,4);
tmp = align_adj{x}([find(time == -1):find(time == 2)],flip(idx(b)))';
h2 = heatmap(tmp);
h2.Title = 'ACh3.0 DMS Example';
h2.XLabel = 'Latency to reward (s)'; h.YLabel = 'Trial';
h2.GridVisible = 'off';
h2.Colormap = parula; h2.ColorLimits = [-6 10];

%% UNIQUE ANALYSIS
% modACh DLS: -3, 7
% allDMS: -1, 7
thres = -2; % CHANGE -- threshold for peaks/pauses to exceed, after baseline
width = 7; % CHANGE -- require that peaks/pauses are at least X samples long (1 sample = 20ms)
thres_used = []; width_used = [];
prc_captured = [];
rew_valMax = cell(length(beh),1); rew_dur = cell(length(beh),1);
imm_valMax = cell(length(beh),1); imm_dur = cell(length(beh),1);
freq = [];
% rew_mu = []; rew_sigma = []; rew_num = [];
% nonrew_freq = []; nonrew_dur = {};
ach2achpause = cell(length(beh),2);

h = waitbar(0, 'pauses');
for x = 1:length(beh); y = 1;
 
    fp = beh(x).FP{1}; Fs = beh(x).Fs; % extract photometry signal from structure
    fp_mu = nanmean(fp); % mean of entire photometry signal
    fp = fp - fp_mu; % SUBTRACT baseline% from fp%, now centered on BASELINE
    rew_true = beh(x).reward(~isnan(align_full{x,y}(1,:))); % excluded non-rewarded trials
    idx_rew = extractEventST([1:length(fp)]', floor(rew_true), floor(rew_true)+100, 1); % identify sample during reward
    idx_mov = extractEventST([1:length(fp)]', beh(x).on, beh(x).off, 1); % identify sample during locomotion
    idx_imm = extractEventST([1:length(fp)]', beh(x).onRest, beh(x).offRest, 1); % identify sample during locomotion
    idx_imm_nonRew = idx_imm(~ismember(idx_imm, idx_rew)); % exclude reward, include rest
    fp_imm = fp(idx_imm_nonRew);
 
    %% Variable threshold and width
    a_rew = align_full{x,y}([find(time==0):find(time==1)],:);
    t_rew = time(find(time==0):find(time==1));
    % a_rew = a_rew - nanmean(align_full{x,y}([find(time == -6):find(time == -2)],:)); % Subtract baseline
    a_rew = a_rew - fp_mu;
    a_rew = a_rew(:,all(~isnan(a_rew))); % remove NaN (unrewarded) trials
    a_rew = a_rew([5:end],:); t_rew = t_rew(5:end); % pause window: [0.2 1]s
    % [m, i] = min(a_rew); % Minima within pause window: [0.2 1]s
    % thres_adj = 0.5*mean(m); % Variable threshold is half max
    
    a = nanmean(a_rew,2);
    thres_adj = 0.5*min(a);
    b = find(a < thres_adj); % Find locations in average reward signal that exceed new threshold
    c = b(1); % First value that crosses threshold
    [~, d] = min(a);
    e = find(a(d:end) > thres_adj); % Find locations that exceed new threshold, after maximum deflection point
    f = e(1) + d-1; % First value that crosses threshold after maximum deflection point
    % width_adj = 0.75*(f-c); % Variable width at half max
    thres_adj = thres;
    width_adj = width;
    % rew_valMax{x} = m(:); % STORE
    
    %%
    prc_captured(x) = 0; counter = 0;
    while prc_captured(x) < 0.9 % Adjust threshold parameters until achieve 90%
        if counter >= 1; thres_adj = thres_adj + 0.1; end
        fp_rew = fp(idx_rew); % photometry during onset window: [-0.2 1]
        [rew_idxMax, rew_idxGood, rew_valMaxAll, rew_crossGroups, rew_crossStart] = findIdxPause(fp_rew, 'pause', thres_adj, width_adj); 
        rew_startGood = rew_crossStart(rew_idxGood,:); % Extract start/end times for good pauses
        rew_valGood = rew_valMaxAll(rew_idxGood);
        tmp = [1:101:length(fp_rew)]'; % Reward start index, in index values of fp_rew
        tmp = ceil(rew_startGood./101); % Divide by length of window to find which reward this pause is in
        tmp((diff(tmp,1,2) ~= 0),:) = nan; % NaN if pause start/end not within same reward window
        uni = unique(tmp); uni = uni(~isnan(uni)); % Unique reward windows that have crossings
        rew_valGoodUnique = nan(length(uni),1); % Initialize
        rew_startGoodUnique = nan(length(uni),2); % Initialize
        rew_idxMaxUnique = nan(length(uni),1); % Initialize
        for c = 1:length(uni)
            cc = find(tmp(:,1) == uni(c)); % idxGood that matches unique reward windows
            [rew_valGoodUnique(c),ii] = min(rew_valGood(cc)); % Keep minimum for multiple crosses in single reward period
            rew_startGoodUnique(c,:) = rew_startGood(cc(ii),:); 
            rew_idxMaxUnique(c) = rew_idxMax(cc(ii));
        end
        rew_valMax{x} = rew_valGoodUnique;
        prc_captured(x) = length(rew_idxMaxUnique)/size(a_rew,2);
        counter = counter + 1; % Add to counter
        if thres_adj > -0.5; thres_adj = -0.5; break; end
    end
    thres_used(x) = thres_adj; width_used(x) = width_adj;
%     pause_tmp = [];
%     for c = 1:length(rew_crossGroups)
%         pause_tmp(c) = min(fp_rew(rew_crossGroups{c})); % Identify minima of threshold crossing groups
%     end 
%     rew_idxMatch = ismember(pause_tmp(:), m(:)); % Where is max onset
    % pause_group = cross_cell(where); 
    % rew_crossStart = rew_crossStart(rew_idxMatch,:); % Extract grouped crossings that match maximum onset group
    
    %%
    [imm_idxMax, imm_idxGood, imm_valMax{x}, ~, imm_crossStart] = findIdxPause(fp_imm, 'pause', thres_adj*2.5, width_adj*0.5); 

    %% ach to pause
    % idxpause_imm = findIdxPause(fp(idx_imm_nonRew),'pause',thres,width);
    % idxpause_rew = findIdxPause(fp(idx_rew),'pause',thres,width);
    [ach2achpause{x,1},sta_time] = getSTA(fp(idx_imm_nonRew), imm_idxMax/Fs, Fs, [-2 2]);
    ach2achpause{x,2} = getSTA(fp(idx_rew), rew_idxMaxUnique/Fs, Fs, [-2 2]);
    
    %%
    rew_dur{x} = []; imm_dur{x} = []; rew_valMax{x} = []; imm_valMax{x} = [];
    for z = 1:2
        if size(ach2achpause{x,z},2) <= 3; continue; end
        a = nanmean(ach2achpause{x,z}(find(sta_time == -0.4):find(sta_time == 0.4),:),2);
        [b,c] = min(a);
        d = find(a < 0.5*b); % Find index at half max
        d = d(1); % First value that crosses threshold
        e = find(a(c:end) > 0.5*b); % Find locations that exceed new threshold, after maximum deflection point
        e = e(1) + c-1; % First value that crosses threshold after maximum deflection point
        if z == 1; imm_dur{x} = e-d; imm_valMax{x} = b;
        elseif z == 2; rew_dur{x} = e-d; rew_valMax{x} = b;
        end
    end
        
%%
    % rew_dur{x} = diff(rew_startGoodUnique,1,2); % Duration of onset peak that exceeds threshold
    % imm_dur{x} = diff(imm_crossStart(imm_idxGood,:),1,2); % Duration of pause/peak that exceeds threshold
    freq(x) = (1/mean(diff(imm_idxMax)))*Fs; % Frequency of maximum deflections
    
    % rew_mu(x,:) = [nanmean(rew_valMax{x}), nanmean(imm_valMax{x})]; % MEAN of pause max amplitude
    % rew_sigma(x,:) = [nanstd(rew_valMax{x}), nanstd(imm_valMax{x})]; % STD of pause max amplitude
    % rew_num(x,:) = [length(find(~isnan(rew_valMax{x}))), length(imm_idxGood)]; % Number of pauses satisfying conditions
    % nonrew_freq(x) = (1/mean(diff(imm_idxMax)))*Fs; % Frequency of maximum deflections
    
    %%
    waitbar(x/length(beh),h);
%%
end
close(h);
 
%% PLOT: Pause Amplitude errorbars
% figure; hold on; 
% errorbar(rew_mu(:,1),rew_sigma(:,1),'c'); errorbar(rew_mu(:,2),rew_sigma(:,2),'k')
% ylabel('Pause Amplitude (%dF/F)'); legend({'rew','non-rew'});
% title(sprintf('Pause Amplitude (thres = %1.1f) (p = %1.2f)',thres,signrank(rew_mu(:,1), rew_mu(:,2))))
 
%% PLOT: HISTOGRAM
% figure;
% for x = 1:length(beh)
%     sp(x) = subplot(4,6,x); hold on
%     %histogram(pause_rew{x},'BinWidth',0.5,'Normalization','probability','FaceColor','b','FaceAlpha',0.2,'EdgeAlpha',0.1);
%     %if ~isempty(pause_nonRew{x})
%     %yyaxis right; histogram(pause_nonRew{x},'BinWidth',0.5,'Normalization','probability','FaceColor','k','FaceAlpha',0.2,'EdgeAlpha',0.1);
%     %end
%     violinplot(nonrew_dur{x});
%     title(sprintf('%s',beh(x).rec));
% end
 
%% n = X mice
tmp = {}; for x = 1:length(beh); tmp{x} = strtok(beh(x).rec,'_'); end
uni = unique(tmp); nAn = length(uni);
mag_an_rew = cell(nAn,1); mag_an_nonrew = cell(nAn,1); % Initialize cell arrays for by animal grouping
dur_an_rew = cell(nAn,1); dur_an_nonrew = cell(nAn,1);
ach2pause_an = cell(2,1);
freq_an = [];
for x = 1:nAn
    ii = find(strcmp(tmp,uni{x}));
    ach2pause_an{1}(:,x) = nanmean([ach2achpause{ii,1}],2);
    ach2pause_an{2}(:,x) = nanmean([ach2achpause{ii,2}],2);
    for y = 1:length(ii)
        mag_an_rew{x} = [mag_an_rew{x}; rew_valMax{ii(y)}];
        mag_an_nonrew{x} = [mag_an_nonrew{x}; imm_valMax{ii(y)}];
        dur_an_rew{x} = [dur_an_rew{x}; rew_dur{ii(y)}];
        dur_an_nonrew{x} = [dur_an_nonrew{x}; imm_dur{ii(y)}];
    end
    freq_an(x) = nanmean(freq(ii));
end

mag_an_mu = [cellfun(@nanmean, mag_an_rew), cellfun(@nanmean, mag_an_nonrew)]; 
dur_an_mu = [cellfun(@nanmean, dur_an_rew), cellfun(@nanmean, dur_an_nonrew)]; 
dur_an_mu = dur_an_mu.*(1/Fs); % convert to seconds

%% PLOT: n = X mice -- Pause MAGNTIDUE + DURATION errorbars
fig = figure; fig.Position(3) = 1375;

subplot(1,3,1); hold on; clr = {'k','b'};
% plot([0 0],[-6 3],'--k'); plot([-1 1],[-2 -2],'--k');  
for y = 1:length(ach2pause_an)
shadederrbar(sta_time, nanmean(ach2pause_an{y},2), SEM(ach2pause_an{y},2), clr{y}); hold on
end
xlabel('Latency to ACh pause (s)'); xlim([-1 1]); xticks([-1:0.5:1]);
ylabel('ACh3.0 (%dF/F)'); ylim([-4 2]); yticks([-6:4]);
axis('square')
title(sprintf('ACh to ACh (n = %d) (every %1.2f s)',nAn,1/nanmean(freq_an)));

subplot(1,3,2); hold on
%
plot_me = mag_an_mu; % CHANGE
%
p = signrank(plot_me(:,1),plot_me(:,2));%,'tail','left'); % null hypothesis: mean difference = 0 VS difference < 0 (reward pause MAGNITUDE is more positive than non reward pause)
errorbar(abs(nanmean(plot_me)'),abs(SEM(plot_me,1)'),'.k','MarkerSize',20);
plot([1.15; 1.85].*ones(2,nAn),abs([plot_me']),'.:k','MarkerSize',20);
ylabel('magnitude (%dF/F)'); ylim([0 4]); yticks([0:4]);
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'rew','immobility'});
axis('square')
title(sprintf('Pause Magnitude (thres = %d) (p = %1.2f)',thres,p));

subplot(1,3,3); hold on
% PLOT: n = X mice -- Pause Duration errorbars
plot_me = dur_an_mu.*1000; % CHANGE
%
p = signrank(plot_me(:,1),plot_me(:,2)); %,'tail','right'); % null hypothesis: mean difference = 0 VS difference < 0 (reward pause DURATION is more positive than non reward pause)
errorbar(abs(nanmean(plot_me)'),abs(SEM(plot_me,1)'),'.k','MarkerSize',20);
plot([1.15; 1.85].*ones(2,nAn),abs([plot_me']),'.:k','MarkerSize',20);
ylabel('duration (ms)'); ylim([0 400]); yticks([0:100:400]);
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'rew','immobility'});
axis('square')
title(sprintf('Pause Duration (width = %d) (p = %1.2f)',width,p));

movegui(gcf,'center');
%% PLOT: n = X mice -- HISTOGRAM
% figure;
% for x = 1:nAn
%     sp(x) = subplot(3,4,x); hold on
%     % histogram(pause_an_rew{x},'BinWidth',0.5,'Normalization','probability','FaceColor','b','FaceAlpha',0.2,'EdgeAlpha',0.1);
%     % yyaxis right; histogram(pause_an_nonrew{x},'BinWidth',0.5,'Normalization','probability','FaceColor','k','FaceAlpha',0.2,'EdgeAlpha',0.1);
%     violinplot(dur_an_nonrew{x});
%     title(sprintf('%s',uni{x}));
% end
% linkaxes(sp,'y');

%% PLOT: ACh to ACh pause
% figure; hold on; clr = {'k','b'};
% plot([0 0],[-6 3],'--k'); plot([-1 1],[-2 -2],'--k');  
% for y = 1:length(ach2pause_an)
% shadederrbar(sta_time, nanmean(ach2pause_an{y},2), SEM(ach2pause_an{y},2), clr{y}); hold on
% end
% xlabel('Latency to ACh pause (s)'); xlim([-1 1]); xticks([-1:0.5:1]);
% ylabel('ACh3.0 (%dF/F)'); yticks([-6:2:3]);
% axis('square')
% title(sprintf('ACh to ACh pause (n = %d mice)',nAn));
