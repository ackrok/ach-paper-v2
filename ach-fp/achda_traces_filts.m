x = 12;

fp_1 = beh(x).FP{y(1)}; Fs = beh(x).Fs; % extract photometry signal from structure
fp_mu = nanmean(fp_1); % mean of entire photometry signal
fp_1 = fp_1 - fp_mu; % SUBTRACT baseline% from fp%, now centered on BASELINE
if isfield(beh,'reward')
    idx_rew = extractEventST([1:length(fp_1)]', floor(beh(x).reward)-100, floor(beh(x).reward)+100, 1); % identify sample during reward
else; idx_rew = []; end
idx_mov = extractEventST([1:length(fp_1)]', beh(x).on, beh(x).off, 1); % identify sample during locomotion
idx_imm = extractEventST([1:length(fp_1)]', beh(x).onRest, beh(x).offRest, 1); % identify sample during locomotion
idx_imm_nonRew = idx_imm(~ismember(idx_imm, idx_rew)); % exclude reward, include rest
% fp_1_imm = fp_1(idx_imm_nonRew); 
fp_2 = beh(x).FP{y(2)}; % dopamine
fp_2 = fp_2 - nanmean(fp_2);
% fp_2_imm = fp_2(idx_imm_nonRew);

fig = figure; 
sp(1) = subplot(2,2,1); hold on; plot(fp_1(idx_imm_nonRew),'g'); plot(fp_2(idx_imm_nonRew),'m'); title(sprintf('%s',beh(x).rec));

f = [0.1 0.5];
    fp_1_filt = filterFP(fp_1,Fs,f,10,'bandpass');
    fp_2_filt = filterFP(fp_2,Fs,f,10,'bandpass');
    fp_1_imm = fp_1_filt(idx_imm_nonRew); 
    fp_2_imm = fp_2_filt(idx_imm_nonRew); 
sp(2) = subplot(2,2,2); hold on; plot(fp_1_imm,'g'); plot(fp_2_imm,'m'); title('0.1-0.5Hz');

f = [0.5 4];
    fp_1_filt = filterFP(fp_1,Fs,f,10,'bandpass');
    fp_2_filt = filterFP(fp_2,Fs,f,10,'bandpass');
    fp_1_imm = fp_1_filt(idx_imm_nonRew); 
    fp_2_imm = fp_2_filt(idx_imm_nonRew); 
sp(3) = subplot(2,2,3); hold on; plot(fp_1_imm,'g'); plot(fp_2_imm,'m'); title('0.5-4Hz');

f = [4 8];
    fp_1_filt = filterFP(fp_1,Fs,f,10,'bandpass');
    fp_2_filt = filterFP(fp_2,Fs,f,10,'bandpass');
    fp_1_imm = fp_1_filt(idx_imm_nonRew); 
    fp_2_imm = fp_2_filt(idx_imm_nonRew); 
sp(4) = subplot(2,2,4); hold on; plot(fp_1_imm,'g'); plot(fp_2_imm,'m'); title('4-8Hz');
linkaxes(sp,'x');
xlim([2e4 2.025e4])

