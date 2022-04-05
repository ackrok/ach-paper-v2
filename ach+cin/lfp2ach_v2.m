[fName, fPath] = uigetfile('*.mat','MultiSelect','On');

%%
mat = struct;
h = waitbar(0, 'tuning unit to LFP');
for x = 1:length(fName)
    load(fullfile(fPath,fName{x}));

    ii = find(strcmp({behACh.rec},lfp.rec));
    if isempty(ii); continue; end
    mat(x).rec = lfp.rec;
    fp = behACh(ii).FP{1}; Fs = 50;
    fp = fp - nanmean(fp);
    idx_imm = extractEventST([1:length(fp)]', behACh(ii).onRest, behACh(ii).offRest, 1); % identify sample during locomotion
    fpImm = fp(idx_imm);
    fpImm = fpImm - nanmean(fpImm);
    
    sig = lfp.data(:);
    sigFs = lfp.samplingRate;
    sigDown = downsampleTLab(sig,sigFs/Fs,2); 
    sigImm = sigDown(idx_imm);
    
    %%
    [mat(x).corr, lags] = xcorr(fpImm, sigImm, 2*Fs, 'coeff');
    tmp_shuff = []; new_sig = sigImm;
    for s = 1:50
        % new_sig = circshift(new_sig, Fs);
        % tmp_shuff(:,s) = xcorr(fp, new_sig, 10*Fs, 'coeff');
        % tmp_shuff(:,s) = xcorr(fp_imm(randperm(length(fp_imm))), new_sig(randperm(length(new_sig))), 10*Fs, 'coeff');
        tmp_shuff(:,s) = xcorr(fpImm, new_sig(randperm(length(new_sig))'), 2*Fs, 'coeff');
    end
    mat(x).corr_5 = prctile(tmp_shuff, 5, 2);
    mat(x).corr_50 = prctile(tmp_shuff, 50, 2);
    mat(x).corr_95 = prctile(tmp_shuff, 95, 2);
    
    %%
    thres_peak = 2; % CHANGE -- threshold for peaks/pauses to exceed, after baseline
    thres_pause = -1;
    width_peak = 5; % CHANGE -- require that peaks/pauses are at least X samples long (1 sample = 20ms)
    locs = findIdxPause(fp,'peak',thres_peak,width_peak);
    locs = extractEventST(locs, behACh(ii).onRest, behACh(ii).offRest, 1);
    locs = locs./Fs;
    [mat(x).staPeak, sta_time] = getSTA(double(sig), locs, sigFs, [-1 1]);
    
    locs = findIdxPause(fp,'pause',thres_pause,width_peak);
    locs = extractEventST(locs, behACh(ii).onRest, behACh(ii).offRest, 1);
    locs = locs./Fs; % convert to seconds
    [mat(x).staPause, ~] = getSTA(double(sig), locs, sigFs, [-1 1]);
    
    waitbar(x/length(fName),h);
end; close(h);

staPeak = []; staPause = [];
for x = 1:length(mat)
    staPeak(:,x) = nanmean(mat(x).staPeak,2);
    staPause(:,x) = nanmean(mat(x).staPause,2);
end

%%
fig = figure; fig.Position(2) = 1000;
subplot(1,2,1); hold on
shadederrbar(lags/Fs, nanmean([mat(x).corr],2), SEM([mat(x).corr],2), 'b');
shadederrbar(lags/Fs, nanmean([mat(x).corr_50],2), SEM([mat(x).corr_5],2), 'k');
xlabel('Lags (s)');
ylabel('Coefficient'); ylim([-0.2 0.1])
title('cross corr: LFP vs ACh'); axis('square');

subplot(1,2,2); hold on
shadederrbar(sta_time, nanmean(staPeak,2), SEM(staPeak,2), 'c');
shadederrbar(sta_time, nanmean(staPause,2), SEM(staPause,2), 'm');
xlabel('Latency (s)'); xlim([-0.3 0.3]);
ylabel('LFP'); 
title('LFP to ACh pause/peak'); axis('square');

%%
figure;
for x = 1:4
    subplot(2,2,x); hold on
    shadederrbar(sta_time, nanmean(mat(x).staPeak,2), SEM(mat(x).staPeak,2), 'c');
    shadederrbar(sta_time, nanmean(mat(x).staPause,2), SEM(mat(x).staPause,2), 'm');
end