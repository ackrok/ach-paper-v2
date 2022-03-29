% x = 2; 
% data.gen.acqFs = 2000; 
% data.acq.nFPchan = 2; data.acq.FP = beh(x).rawFP; data.acq.FPnames = beh(x).FPnames;
% data.acq.wheel = beh(x).wheel; 
% data.acq.time = [1/2000:1/2000:length(beh(x).wheel)/2000]';
% 
% data = processFP(data, params);
% data = processMov(data, params);
% data = processOnsetOffset(data, params);
% data = processRestOnsetOffset(data, params);
% 
% beh(x).FP = data.final.FP; beh(x).Fs = 50;
% beh(x).vel = data.final.vel;
% beh(x).on = data.final.mov.onsets; beh(x).off = data.final.mov.offsets;
% beh(x).onRest = data.final.mov.onsetsRest; beh(x).offRest = data.final.mov.offsetsRest;
load('R:\homes\ack466\ACh paper\6_ACh+DA\gcamp_str+snc.mat');

x = 2; 
figure; hold on
clr = {'b','r'};
for y = 1:2; plot(beh(x).time, beh(x).FP{y} - nanmean(beh(x).FP{y}), clr{y}); end
plot(beh(x).time, getAcc(beh(x).vel), 'k');
legend({'STR','SNc','vel'}); 
title(sprintf('%s',beh(x).rec)); ylabel('Fluorescence (%dF/F)')
xlim([500 525]); xlabel('Time(s)')

%%
corr_mat = nan(1001,length(beh));
corr_5 = nan(1001,length(beh)); corr_50 = corr_5; corr_95 = corr_5; 

for x = 1:length(beh)
    fp_1 = beh(x).FP{1}; Fs = beh(x).Fs; % extract photometry signal from structure
    fp_mu = nanmean(fp_1); % mean of entire photometry signal
    fp_1 = fp_1 - fp_mu; % SUBTRACT baseline% from fp%, now centered on BASELINE
    idx_mov = extractEventST([1:length(fp_1)]', beh(x).on, beh(x).off, 1); % identify sample during locomotion
    idx_imm = extractEventST([1:length(fp_1)]', beh(x).onRest, beh(x).offRest, 1); % identify sample during locomotion
    fp_1_imm = fp_1(idx_imm); 
    
    fp_2 = beh(x).FP{2}; % dopamine
    fp_2 = fp_2 - nanmean(fp_2);
    fp_2_imm = fp_2(idx_imm);
    
    
    %%
    if isempty(fp_1_imm) || isempty(fp_2_imm); continue; end

    %% cross-correlation, whole signal
    [corr_mat(:,x),lags] = xcorr(fp_1_imm, fp_2_imm, 10*Fs, 'coeff');
    tmp_shuff = []; new_fp_2_imm = fp_2_imm;
    for s = 1:50
        new_fp_2_imm = circshift(new_fp_2_imm, Fs);
        tmp_shuff(:,s) = xcorr(fp_1_imm, new_fp_2_imm, 10*Fs, 'coeff');
    end
    corr_5(:,x) = prctile(tmp_shuff, 5, 2);
    corr_50(:,x) = prctile(tmp_shuff, 50, 2);
    corr_95(:,x) = prctile(tmp_shuff, 95, 2);
end

%% 
corr_adj = corr_mat - nanmean(corr_mat([1:50],:));
corr_adj = corr_mat;

fig = figure; fig.Position(3) = 1000;
for x = 1:2 
    sp(x) = subplot(1,2,x); hold on
    plot([0 0],[-0.1 0.5],'--k');
    plot(lags/Fs, corr_adj(:,x), 'b');
    shadederrbar(lags/Fs, nanmean(corr_50(:,x))*ones(length(lags),1), nanmean(corr_5(:,x))*ones(length(lags),1), 'k');
    xlabel('Lag (s)'); xlim([-5 5]);
    ylabel('Sample Cross Correlation');
    title(sprintf('%s',beh(x).rec));
    % Interpretation:
    % peak corr has POSITIVE lag, so DA signal is equal to ACh signal shifted
    % by 10 samples (200ms) to the left
end