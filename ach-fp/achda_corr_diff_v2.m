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

beh = b2flox;

beh = allDMS; beh(13) = [];

%%
mat = struct;

h = waitbar(0, 'cross corr');
for x = 1:length(beh); y = [1 2]; %CHANGE - which FP signal to run CCG or ACG on
    
    %% acetylcholine
    fp_1 = beh(x).FP{y(1)}; Fs = beh(x).Fs; % extract photometry signal from structure
    fp_mu = nanmean(fp_1); % mean of entire photometry signal
    fp_1 = fp_1 - fp_mu; % SUBTRACT baseline% from fp%, now centered on BASELINE
    fp_2 = beh(x).FP{y(2)}; % dopamine
    fp_2 = fp_2 - nanmean(fp_2);
    
    if isfield(beh,'reward')
        idx_rew = extractEventST([1:length(fp_1)]', floor(beh(x).reward)-100, floor(beh(x).reward)+100, 1); % identify sample during reward
    else; idx_rew = []; end
    idx_mov = extractEventST([1:length(fp_1)]', beh(x).on, beh(x).off, 1); % identify sample during locomotion
    idx_mov_nonRew = idx_mov(~ismember(idx_mov, idx_rew)); % exclude reward, include locomotion
    idx_imm = extractEventST([1:length(fp_1)]', beh(x).onRest, beh(x).offRest, 1); % identify sample during locomotion
    idx_imm_nonRew = idx_imm(~ismember(idx_imm, idx_rew)); % exclude reward, include rest
    if isempty(idx_imm_nonRew); continue; end
    
    %% bandpass filter
    f = [0.5 4];
    fp_1_filt = filterFP(fp_1,Fs,f,10,'bandpass');
    fp_2_filt = filterFP(fp_2,Fs,f,10,'bandpass');
    
    %% cross-correlation IMMOBILITY
    fp_1_sub = fp_1_filt(idx_imm_nonRew); fp_1_sub = fp_1_sub - nanmean(fp_1_sub);
    fp_2_sub = fp_2_filt(idx_imm_nonRew); 
    fp_2_sub = [fp_2_sub(1); diff(fp_2_sub)]; % first derivative
    
    % [xcf, lags, bounds] = crosscorr(fp_ach_imm, fp_da_imm,'NumLags',100,'NumSTD',3);
    % corr_achda(:,x) = xcf; corr_achda_bounds(:,x) = bounds; % STORE
    % [corr_achda_shuff(:,x),~,~] = crosscorr(fp_ach_imm(randperm(length(fp_ach_imm))), fp_da_imm(randperm(length(fp_ach_imm))),'NumLags',100,'NumSTD',3);
    [corr_mat,lags] = xcorr(fp_1_sub, fp_2_sub, 10*Fs, 'coeff');
    tmp_shuff = []; new_fp_2_sub = fp_2_sub;
    for s = 1:50
        new_fp_2_sub = circshift(new_fp_2_sub, Fs);
        tmp_shuff(:,s) = xcorr(fp_1_sub, new_fp_2_sub, 10*Fs, 'coeff');
        % tmp_shuff(:,s) = xcorr(fp_1_imm(randperm(length(fp_1_imm))), fp_2_imm(randperm(length(fp_2_imm))), 10*Fs, 'coeff');
        % tmp_shuff(:,s) = xcorr(fp_1_imm, fp_2_imm(randperm(length(fp_2_imm))), 10*Fs, 'coeff');
    end
    mat(x).rec = beh(x).rec;
    mat(x).corr_imm = corr_mat;
    mat(x).shuff_imm = prctile(tmp_shuff, [5 50 95], 2);
    
    %% cross-correlation LOCOMOTION
    fp_1_sub = fp_1_filt(idx_mov_nonRew); fp_1_sub = fp_1_sub - nanmean(fp_1_sub);
    fp_2_sub = fp_2_filt(idx_mov_nonRew); 
    fp_2_sub = [fp_2_sub(1); diff(fp_2_sub)]; % first derivative
    
    [corr_mat,~] = xcorr(fp_1_sub, fp_2_sub, 10*Fs, 'coeff');
    tmp_shuff = []; new_fp_2_sub = fp_2_sub;
    for s = 1:50
        new_fp_2_sub = circshift(new_fp_2_sub, Fs);
        tmp_shuff(:,s) = xcorr(fp_1_sub, new_fp_2_sub, 10*Fs, 'coeff');
    end
    mat(x).corr_mov = corr_mat;
    mat(x).shuff_mov = prctile(tmp_shuff, [5 50 95], 2);
    
    %% cross-correlation REWARD
    if ~isempty(idx_rew)
        fp_1_sub = fp_1_filt(idx_rew); fp_1_sub = fp_1_sub - nanmean(fp_1_sub);
        fp_2_sub = fp_2_filt(idx_rew); 
        fp_2_sub = [fp_2_sub(1); diff(fp_2_sub)]; % first derivative
        
        [corr_mat,~] = xcorr(fp_1_sub, fp_2_sub, 10*Fs, 'coeff');
        tmp_shuff = []; new_fp_2_sub = fp_2_sub;
        for s = 1:50
            new_fp_2_sub = circshift(new_fp_2_sub, Fs);
            tmp_shuff(:,s) = xcorr(fp_1_sub, new_fp_2_sub, 10*Fs, 'coeff');
        end
        mat(x).corr_rew = corr_mat;
        mat(x).shuff_rew = prctile(tmp_shuff, [5 50 95], 2);
        mat(x).a = 1;
    else; mat(x).a = 0;
    end
    
%%
    waitbar(x/length(beh),h);
end
close(h);

%% N = X mice
rec = {}; for x = 1:length(beh); rec{x} = strtok(beh(x).rec,'-'); end
uni = unique(rec); nAn = length(uni);
corr_an = cell(3,1); corr_an_5 = cell(3,1); corr_an_50 = cell(3,1);
for x = 1:nAn
    ii = find(strcmp(rec,uni{x})); % index of matching recordings for this animal
    y = 1;
    tmp = [mat(ii).corr_imm];
    corr_an{y} = [corr_an{y}, nanmean(tmp,2)];
    tmp = [mat(ii).shuff_imm];
    corr_an_5{y} = [corr_an_5{y}, nanmean(tmp(:,1:3:size(tmp,2)),2)];
    corr_an_50{y} = [corr_an_50{y}, nanmean(tmp(:,2:3:size(tmp,2)),2)];
    
    y = 2;
    tmp = [mat(ii).corr_mov]; 
    corr_an{y} = [corr_an{y}, nanmean(tmp,2)];
    tmp = [mat(ii).shuff_mov];
    corr_an_5{y} = [corr_an_5{y}, nanmean(tmp(:,1:3:size(tmp,2)),2)];
    corr_an_50{y} = [corr_an_50{y}, nanmean(tmp(:,2:3:size(tmp,2)),2)];
    
    y = 3;
    ii_rew = ii(ismember(ii, find([mat.a]))); % retain index of recordings that include reward
    tmp = [mat(ii_rew).corr_rew]; 
    corr_an{y} = [corr_an{y}, nanmean(tmp,2)];
    tmp = [mat(ii_rew).shuff_rew];
    corr_an_5{y} = [corr_an_5{y}, nanmean(tmp(:,1:3:size(tmp,2)),2)];
    corr_an_50{y} = [corr_an_50{y}, nanmean(tmp(:,2:3:size(tmp,2)),2)];
end

%% PLOT AVG
fig = figure; fig.Position(3) = 1375;
clr = {'r','g','b'}; lbl = {'IMM','MOV','REW'};
for y = 1:3
sp(y) = subplot(1,3,y); hold on
plot([0 0],[-0.4 0.8],'--k');
shadederrbar(lags/Fs, nanmean(corr_an_50{y},2), nanmean(corr_an_5{y},2), 'k'); hold on
plot(lags/Fs, corr_an{y}, 'Color', [0 0 0 0.1]);
shadederrbar(lags/Fs, nanmean(corr_an{y},2), SEM(corr_an{y},2), clr{y});
xlabel('Lag (s)'); xlim([-1.5 1.5]);
ylabel('Sample Cross Correlation'); ylim([-0.5 1]); yticks([-0.5:0.5:1]);
title(sprintf('ACh/diff(DA) (n = %d mice) %s',nAn,lbl{y})); axis('square');
end
linkaxes(sp,'y'); linkaxes(sp,'x');
movegui(gcf,'center');
% Interpretation:
% peak corr has POSITIVE lag, so DA signal is equal to ACh signal shifted
% by 10 samples (200ms) to the left
%
% if peak corr has negative lag, then means that xcorr(vec1, vec2)
% = vec2 is equal to vec1 shifted by X samples to the right 
% = vec2 follows vec1 by X samples

%% HEATMAP
win = [-1 1];
fig = figure; fig.Position(3) = 1375;
lbl = {'IMM','MOV','REW'};
for y = 1:3
sp(y) = subplot(1,3,y);
a = corr_an{y};
b = a(lags == 0,:); % find delta @ lag = 0
%[c, ii] = sort(b); % sort in ascending order
h(y) = heatmap(a([find(lags == win(1)*Fs) : find(lags == win(2)*Fs)],:)');
h(y).Title = sprintf('ACh/diff(DA) %s',lbl{y});
h(y).XLabel = 'Lags (s)'; % h(y).YLabel = 'Animal';
h(y).Colormap = jet; h(y).GridVisible = 'off';
h(y).ColorLimits = [-0.5 1];
end
movegui(gcf,'center');

%% COMPARE STATS
corr_mat = []; max_val = []; max_lag = [];
for y = 1:3
    corr_mat(:,y) = nanmean(corr_an{y},2);
    [max_val(:,y),ii] = max(corr_an{y}); % maximum correlation coefficient
    max_lag(:,y) = lags(ii)/Fs*1000; % lag of maximum, in milliseconds
end

fig = figure; fig.Position(3) = 1375;
subplot(1,3,1); hold on
for y = 1:3; plot(lags/Fs, corr_mat(:,y), clr{y}); end
plot([0 0],[-0.4 0.8],'--k','DisplayName','t=0');
xlabel('Lag (s)'); xlim([-1.5 1.5]);
ylabel('Correlation Coefficient'); ylim([-0.5 1]); legend(lbl);
title(sprintf('ACh/diff(DA) (n = %d mice)',nAn)); axis('square');

subplot(1,3,2); hold on
a = max_val;
plot(a','--','Color',[0 0 0 0.5]);
errorbar(nanmean(a),SEM(a,1),'.-r','MarkerSize',20);
xlim([0.5 3.5]); xticks([1:3]); xticklabels(lbl);
ylabel('MAX corr coeff'); ylim([0 1]); yticks([0:0.2:1]);
title(sprintf('MAX signrank (i/m: %1.3f) (i/r: %1.3f)',signrank(a(:,1),a(:,2)),signrank(a(:,1),a(:,3)))); axis('square');

subplot(1,3,3); hold on
a = max_lag;
plot(a','--','Color',[0 0 0 0.5]);
errorbar(nanmean(a),SEM(a,1),'.-r','MarkerSize',20);
xlim([0.5 3.5]); xticks([1:3]); xticklabels(lbl);
ylabel('LAG max corr coeff (ms)'); ylim([-100 100]); yticks([-100:50:100]);
title(sprintf('LAG signrank (i/m: %1.3f) (i/r: %1.3f)',signrank(a(:,1),a(:,2)),signrank(a(:,1),a(:,3)))); axis('square');

movegui(gcf,'center');