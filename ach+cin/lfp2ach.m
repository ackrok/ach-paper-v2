x = 1;
lfp_band = bandpass(roi(x).lfp(1:length(roi(x).fp)), [0.5 4], roi(x).Fs);

%% (1) LFP band + ACh cross-correlation
corr_mat = []; corr_5 = []; corr_50 = []; corr_95 = [];

a = [2 14];
figure; hold on
for b = 1:length(a); x = a(b);
    if isempty(roi(x).lfp); continue; end
    fp = roi(x).fp; Fs = roi(x).Fs;
    lfp_band = bandpass(roi(x).lfp(1:length(fp)), [0.5 4], Fs);
    [corr_tmp,lags] = xcorr(fp, lfp_band, 10*Fs, 'coeff');
    corr_mat = [corr_mat, corr_tmp];
    tmp_shuff = []; new_lfp = lfp_band;
    for s = 1:50
        new_lfp = circshift(new_lfp, Fs);
        % new_lfp = lfp_band(randperm(length(lfp_band)));
        tmp_shuff(:,s) = xcorr(fp, new_lfp, 10*Fs, 'coeff');
    end
    corr_5  = [corr_5,  prctile(tmp_shuff, 5, 2)];
    corr_50 = [corr_50, prctile(tmp_shuff, 50, 2)];
    corr_95 = [corr_95, prctile(tmp_shuff, 95, 2)];
    
    sp(b) = subplot(1,2,b); hold on
    shadederrbar(lags/Fs, corr_50(:,b), corr_5(:,b), 'k');
    plot(lags/Fs, corr_tmp, 'b');
    title(sprintf('%s: corr LFPband vs ACh',roi(x).rec),'Interpreter','none'); axis('square');
    xlim([-2 2]); xlabel('Lag (s)');
    ylim([-0.1 0.1]); ylabel('corr coeff');
end

% figure; hold on
% plot([0 0],[-0.7 0.2],'--k');
% plot(lags/Fs, corr_ma 'Color', [0 0 0 0.1]);
% shadederrbar(lags/Fs, nanmean(corr_mat,2), SEM(corr_mat,2), 'b');
% xlabel('Lag (s)'); ylabel('Sample Cross Correlation');
% title(sprintf('ACh/DA cross-correlation (n = %d recs)',size(corr_mat,2)));

% figure; 
% for x = 1:size(corr_mat,2)
%     sp(x) = subplot(2,1,x); hold on
%     shadederrbar(lags/Fs, corr_50(:,x), corr_5(:,x), 'k');
%     plot(lags/Fs, corr_mat(:,x), 'b');
% end
linkaxes(sp,'x'); linkaxes(sp,'y');
xlim([-2 2]); ylim([-0.1 0.1]);

%% (2) LFP band aligned to ACh peaks, pauses

for x = [2 14]
    fp = roi(x).fp; Fs = roi(x).Fs;
    fp = fp - nanmean(fp); % subtract mean of baseline recording to center around 0%
    thres_peak = 2; % CHANGE -- threshold for peaks/pauses to exceed, after baseline
    thres_pause = -1;
    width_peak = 5; % CHANGE -- require that peaks/pauses are at least X samples long (1 sample = 20ms)
    locs_peak = findIdxPause(fp,'peak',thres_peak,width_peak);
    locs_pause = findIdxPause(fp,'pause',thres_pause,width_peak);
    
    idx_imm = extractEventST([1:length(fp)]',roi(x).onRest,roi(x).offRest,1);
    locs_peak_imm = locs_peak(ismember(locs_peak, idx_imm));
    locs_pause_imm = locs_pause(ismember(locs_pause, idx_imm));
    
    lfp_band = bandpass(roi(x).lfp(1:length(fp)), [0.5 4], Fs);
    [sta_peak, time] = getSTA(lfp_band, locs_peak_imm./Fs, Fs, [-1 1]);
    sta_pause = getSTA(lfp_band, locs_pause_imm./Fs, Fs, [-1 1]);
    
    roi(x).sta_peak = sta_peak; roi(x).sta_pause = sta_pause;
end

a = [2 14];
fig = figure; fig.Position(3) = 1000; hold on
for b = 1:length(a); x = a(b);
    sp(b) = subplot(1,2,b); hold on
    shadederrbar(time, nanmean(roi(x).sta_peak,2), SEM(roi(x).sta_peak,2), 'c');
    shadederrbar(time, nanmean(roi(x).sta_pause,2), SEM(roi(x).sta_pause,2), 'r');
    plot([0 0],[-100 100],'--k');
    legend({'ACh peak','','ACh pause'});
    ylabel('LFP [0.5 4]Hz band'); ylim([-100 100]);
    xlabel('Latency to ACh pk/ps (s)'); xticks([-1:0.25:1]);
    title(sprintf('%s: LFPband to ACh',roi(x).rec),'Interpreter','none'); axis('square');
end
linkaxes(sp,'x'); linkaxes(sp,'y');

%% LFP band to CIN spikes
a = [2 14];
fig = figure; fig.Position(3) = 1000; hold on
for b = 1:length(a); x = a(b);
    idx = find(strcmp({cinACh.rec},roi(x).rec));
    st = {cinACh(idx).st};
    lfp2cin = cell(length(st),1);
    st_imm = cell(length(st),1);
    
    sp(b) = subplot(1,2,b); hold on
    for y = 1:length(st)
        st_imm{y} = extractEventST(st{y},roi(x).onRest,roi(x).offRest,1);
        [sta_mat, time] = getSTA(roi(x).lfp, st_imm{y}, roi(x).Fs, [-1 1]);
        lfp2cin{y} = sta_mat;
        plot(time, nanmean(lfp2cin{y},2));
    end
    ylabel('LFP [0.5 4]Hz band'); ylim([-100 100]);
    xlabel('Latency to CIN spike (s)'); xticks([-1:0.25:1]);
    title(sprintf('%s: LFPband to CIN st',roi(x).rec),'Interpreter','none'); axis('square');
end

%% EXTRACT LFP
for x = 2
    nChan = 128;
    lfpChan = roi(x).maxChan(1);
    acqFs = 30000;
    dsFs = roi(x).Fs; % match sampling rate to photometry for cross-correlation
    fPath = 'C:\Users\Anya\Desktop\IV_LOCAL\ACh';
    fName = roi(x).rec;

    fprintf('%s loading LFP ... ',fName);
    lfpVec = loadLFP(nChan,lfpChan,acqFs,dsFs); % Load LFP
    fprintf('done ... ');
    roi(x).lfp = lfpVec(:);

    lfp = struct; 
    lfp.mouserec = roi(x).rec;
    lfp.chan = num2str(lfpChan); 
    lfp.Fs = dsFs; 
    lfp.lfp_ds = lfpVec;
    save(['C:\Users\Anya\Desktop\IV_LOCAL\lfp_local\',lfp.mouserec,'_ch',lfp.chan,'_LFP'],'lfp');
    fprintf(' saved \n');
end