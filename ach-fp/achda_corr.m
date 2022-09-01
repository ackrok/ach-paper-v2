%Cross-correlation between ACh and DA photometry signals acquired during
%dual recordings
%
%   Description: This script is for running cross-correlation analysis on
%   two continuous photometry signals (ACh, DA) after isolating time
%   periods when animal is immobile (no movement, no reward)
%
%   Author: Anya Krok, December 2021

%% 
% good = [1:length(modAChDA)]; rmv = [3 5 18 22:24 27:28 30 33:34 36 41:42]; good(rmv) = []; % acceleration
% good = [10:16,20,22:29,32:35,40:41,44:47]; % reward
good = [1:length(modAChDA)]; rmv = 42; good(rmv) = [];
beh = modAChDA(good);

beh = dual;

%%
mat = struct;
corr_cell = cell(3,4);
for a = 1:3; for b = 1:4; corr_cell{a,b} = nan(501,length(beh)); end; end

h = waitbar(0, 'cross corr');
for x = 1:length(beh); y = [1 2]; %CHANGE - which FP signal to run CCG or ACG on
    
    %% extract signals
    fp_mat = [];
    fp_mat(:,1) = beh(x).FP{y(1)}; Fs = beh(x).Fs; % extract photometry signal from structure
    fp_mat(:,1) = fp_mat(:,1) - nanmean(fp_mat(:,1)); % subtract baseline (mean of entire photometry signal) from fp
    fp_mat(:,2) = beh(x).FP{y(2)}; % dopamine
    fp_mat(:,2) = fp_mat(:,2) - nanmean(fp_mat(:,2)); % subtract baseline (mean of entire photometry signal) from fp
    % fp_mat(:,2) = filterFP(fp_mat(:,2),Fs,[0.5 4],10,'bandpass'); 
    
%     fp_filt = [];
%     passband = [0.5 4]; % bandpass filter
%     for y = 2
%         fp_filt(:,y) = filterFP(fp_mat(:,y),Fs,passband,10,'bandpass'); 
%     end
%     testdiff = diff(fp_filt(:,2)); testdiff = testdiff - nanmean(testdiff);
%     fp_filt(:,2) = testdiff;
%     testdiff = cumsum(fp_filt(:,1)); testdiff - nanmean(testdiff);
%     fp_filt(:,1) = testdiff

    %% behavioral states
    if isfield(beh,'reward')
        idx_rew = extractEventST([1:length(fp_mat)]', floor(beh(x).reward), floor(beh(x).reward)+50, 1); % identify sample during reward
    else; idx_rew = []; end
    idx_mov = extractEventST([1:length(fp_mat)]', beh(x).on, beh(x).off, 1); % identify sample during locomotion
    idx_mov_nonRew = idx_mov(~ismember(idx_mov, idx_rew)); % exclude reward, include locomotion
    idx_imm = extractEventST([1:length(fp_mat)]', beh(x).onRest, beh(x).offRest, 1); % identify sample during locomotion
    idx_imm_nonRew = idx_imm(~ismember(idx_imm, idx_rew)); % exclude reward, include rest
    idx_cell = cell(3,1); idx_cell{1} = idx_imm_nonRew; idx_cell{2} = idx_mov_nonRew; idx_cell{3} = idx_rew; % index into cell array for ease of iteration
    
    %%
    mat(x).rec = beh(x).rec;
    
    fp_sub = [];
    for z = 1:3
        if length(idx_cell{z})< 2; continue; end
        fp_sub = fp_mat(idx_cell{z},:); % signal
        % fp_sub = normalize(fp_mat,1,'range'); % normalize [0 1]
        % fp_sub = fp_filt(idx_cell{z},:); % bandpass filtered signal
        
        % [xcf, lags, bounds] = crosscorr(fp_sub(:,1), fp_sub(:,2),'NumLags',100,'NumSTD',3);
        % [shuff,~,~] = crosscorr(fp_sub(randperm(size(fp_sub,1)),1), fp_sub(randperm(size(fp_sub,2)),2),'NumLags',100,'NumSTD',3);
        [corr_tmp, lags] = xcorr(fp_sub(:,1), fp_sub(:,2), 5*Fs, 'coeff'); % cross-correlation
        
        fp_sub_new = fp_sub(:,2);
        tmp_shuff = []; 
        for s = 1:50
            fp_sub_new = circshift(fp_sub_new, Fs);
            % tmp_shuff(:,s) = xcorr(fp_sub(randperm(size(fp_sub,1)),1), fp_sub(randperm(size(fp_sub,2)),2), 10*Fs, 'coeff');
            % tmp_shuff(:,s) = xcorr(fp_sub(:,1), fp_sub(randperm(size(fp_sub,2)),2), 10*Fs, 'coeff');
            tmp_shuff(:,s) = xcorr(fp_sub(:,1), fp_sub_new, 5*Fs, 'coeff');
        end
        
        mat(x).corr{z} = corr_tmp;
        mat(x).shuff{z} = prctile(tmp_shuff, [5 50 95], 2);
        corr_cell{z,1}(:,x) = corr_tmp;       % cross-correlation
        corr_cell{z,2}(:,x) = prctile(tmp_shuff, 5, 2); % shuffle 5th percentile
        corr_cell{z,3}(:,x) = prctile(tmp_shuff, 50, 2); % shuffle 50th percentile
        corr_cell{z,4}(:,x) = prctile(tmp_shuff, 95, 2); % shuffle 95th percentile
    end
    
%%
    waitbar(x/length(beh),h);
end
close(h);

%N = X mice
corr_an = cell(3,4); min_val = []; min_lag = [];
tmp = {}; for x = 1:length(beh); tmp{x} = strtok(beh(x).rec,'-'); end
uni = unique(tmp); nAn = length(uni);

for x = 1:nAn
    idx = strcmp(tmp,uni{x});
    for a = 1:3; for b = 1:4
        corr_adj = corr_cell{a,b};
        if b == 1; corr_adj = corr_adj - nanmean(corr_adj([1:find(lags == -2)],:)); end
        corr_an{a,b}(:,x) = nanmean(corr_adj(:,idx),2);
        end; end
end
for z = 1:3
    [min_val(:,z), ii] = min(corr_an{z,1});
    min_lag(:,z) = lags(ii)./50;
end

%% PLOT: average full cross-correlation between ACh, DA signals
z = 1;

figure; hold on
plot([0 0],[-0.7 0.2],'--k');
%shadederrbar(lags/Fs, nanmean(corr_shuff,2), SEM(corr_shuff,2), 'k'); hold on
plot(lags/Fs, corr_cell{z,1}, 'Color', [0 0 0 0.1]);
shadederrbar(lags/Fs, nanmean(corr_cell{z,1},2), SEM(corr_cell{z,1},2), 'b');
xlabel('Lag (s)'); ylabel('Sample Cross Correlation');
title(sprintf('ACh/DA cross-correlation (n = %d recs)',size(corr_cell{z,1},2)));
% Interpretation:
% peak corr has POSITIVE lag, so DA signal is equal to ACh signal shifted
% by 10 samples (200ms) to the left

%% HEATMAP
z = 1;
win = [-1 1];
figure;
a = corr_an{z,1};
b = a(lags == 0,:); % find delta @ lag = 0
[c, ii] = sort(b); % sort in ascending order
h = heatmap(a([find(lags == win(1)*Fs) : find(lags == win(2)*Fs)],ii)');
h.Title = 'ACh/DA corr IMM';
h.XLabel = 'Lags (s)'; 
h.YLabel = 'Recording';
h.Colormap = jet; h.GridVisible = 'off';
h.ColorLimits = [-1 0.5];

%% PLOT N = X mice
fig = figure; fig.Position(3) = 1375; clr = {'r','g','b'};
for z = 1:3
    subplot(1,3,z); hold on
    plot([0 0],[-0.6 0.6],'--k');
    a = nanmean(corr_an{z,3},2);
    a = a - nanmean(a(find(lags/Fs == -5):find(lags/Fs == -1),:));
    shadederrbar(lags/Fs, a, nanmean(corr_an{z,4}-corr_an{z,3},2), 'k'); hold on
    a = corr_an{z,1}; a = a - nanmean(a(find(lags/Fs == -5):find(lags/Fs == -1),:));
    plot(lags/Fs, a, 'Color', [0 0 0 0.1]);
    shadederrbar(lags/Fs, nanmean(a,2), SEM(corr_an{z,1},2), clr{z});
    xlabel('Lag (s)'); xlim([-1 1]);
    ylabel('Correlation Coefficient'); ylim([-1 0.5]); yticks([-1:0.25:1]);
    title(sprintf('ACh/DA corr (n = %d mice)',nAn)); axis('square');
end
movegui(gcf,'center');

%% PLOT N = X mice + STATS
fig = figure; fig.Position(3) = 1375; 
subplot(1,3,1); hold on; clr = {'r','g','b'};
plot([0 0],[-1 0.6],'--k');
shadederrbar(lags/Fs, nanmean(corr_an{1,3},2), nanmean(corr_an{1,2},2), 'k'); hold on
for z = 1:3
    shadederrbar(lags/Fs, nanmean(corr_an{z,1},2), SEM(corr_an{z,1},2), clr{z});
end
xlabel('Lag (s)'); xlim([-1 1]);
ylabel('Correlation Coefficient'); ylim([-1 0.5]); yticks([-1:0.25:1]);
title(sprintf('ACh/DA corr (n = %d mice)',nAn)); axis('square');

subplot(1,3,2); hold on; clr = {'r','g','b'};
a = min_val;
violinplot(a);
errorbar(nanmean(a,1), SEM(a,1),'.k', 'MarkerSize', 20);
xlim([0.5 3.5]); xticks([1:3]); xticklabels({'imm','mov','rew'}); 
ylabel('Correlation Coefficient'); ylim([-1 0]); yticks([-1:0.2:1]);
title(sprintf('coeff: signrank(i/m: %1.3f | i/r: %1.3f)', signrank(a(:,1),a(:,2)), signrank(a(:,1),a(:,3)))); axis('square');

subplot(1,3,3); hold on; clr = {'r','g','b'};
a = min_lag.*1000;
violinplot(a);
errorbar(nanmean(a,1), SEM(a,1),'.k', 'MarkerSize', 20);
xlim([0.5 3.5]); xticks([1:3]); xticklabels({'imm','mov','rew'}); 
ylabel('Latency to minimum (ms)'); ylim([-250 0]); yticks([-250:50:0]);
title(sprintf('lag: signrank(i/m: %1.3f | i/r: %1.3f)', signrank(a(:,1),a(:,2)), signrank(a(:,1),a(:,3)))); axis('square');
movegui(gcf,'center');

% Interpretation:
% peak corr has POSITIVE lag, so DA signal is equal to ACh signal shifted
% by 10 samples (200ms) to the left
%
% if peak corr has negative lag, then means that xcorr(vec1, vec2)
% = vec2 is equal to vec1 shifted by X samples to the right 
% = vec2 follows vec1 by X samples

%%
figure; 
subplot(1,2,1); violinplot(min_val); 
ylabel('max corr coeff'); % ylim([-1 0]); yticks([-1:0.2:0]);
title(sprintf('ACh/DA coeff: %1.2f',nanmean(min_val)));
subplot(1,2,2); violinplot(min_time*1000); 
ylabel('max corr lag (ms)'); % ylim([0 250]); yticks([0:50:250]);
title(sprintf('ACh/DA lag: %1.1f ms',nanmean(min_time*1000)));

%% PLOT proportion of CCG > 95% CI
above95 = []; below5 = [];
for x = 1:length(uni)
    a = corr_an(:,x);
    above95 = [above95, a > corr_an_95(:,x)]; %binary vector where CCG passed 95% confidence interval
    below5 = [below5, a < corr_an_5(:,x)]; %binary vector where CCG below 5% confidence interval
end

figure; hold on
bar(lags/Fs, 100*sum(above95,2)/size(above95,2),'FaceColor','b','FaceAlpha',0.5);
bar(lags/Fs, -100*sum(below5,2)/size(below5,2),'FaceColor','b','FaceAlpha',0.5);
xlabel('Lag (s)'); ylabel('Prop of CCG > 95% CI (%)'); ylim([-100 100]);
title(sprintf('Proportion of CCG > 95p CI (n = %d mice)',size(above95,2)))

%% PLOT BEHAVIORAL STATES
