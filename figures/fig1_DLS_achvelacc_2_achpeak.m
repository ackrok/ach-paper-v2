a = [1:length(modACh)]; rmv = [3 17 20 22 23 27 30 35:39 41:43 48 50]; a(rmv) = [];
beh = modACh(a); % Extract recordings with reward
[align_full, time] = plot_fp2event(beh,[-4 2],0); % Align photometry to events

%% Make rest threshold more stringent
for x = 1:length(beh)
    vel = beh(x).vel; vel = abs(vel);

    velThres = 0.05;
    Fs = 50;
    minRestTime = params.mov.minRestTime_rest*Fs;
    timeThres = params.mov.timeThres_rest *Fs;
    timeShift = params.mov.timeShift_rest*Fs;
    minRunTime = params.mov.minRunTime_rest * Fs;
    [onsetInd,offsetInd] = getOnsetOffset(-vel,-velThres,minRunTime,minRestTime,1);
    [onsetInd,offsetInd] = adjOnsetOffset(onsetInd,offsetInd,timeThres,vel);
    onsetInd = onsetInd+timeShift; offsetInd = offsetInd-timeShift;
    beh(x).onRest = onsetInd;
    beh(x).offRest = offsetInd;
end

%%
ach2achpeak = cell(length(beh),2); 
vel2achpeak = ach2achpeak; acc2achpeak = vel2achpeak;
thres_keep = [];

h = waitbar(0, 'peaks - onset vs non');
for x = 1:length(beh); y = 1;
    %% 
    thres = 4; % from unique_onset_v2
    width = 6.5; % CHANGE -- require that peaks/pauses are at least X samples long (1 sample = 20ms)
    
    %%
    fp = beh(x).FP{1}; Fs = beh(x).Fs; % extract photometry signal from structure
    fp_mu = nanmean(fp); % mean of entire photometry signal
    fp = fp - fp_mu; % SUBTRACT baseline% from fp%, now centered on BASELINE
    idx_mov = extractEventST([1:length(fp)]', beh(x).on, beh(x).off, 1); % identify sample during locomotion
    idx_imm = extractEventST([1:length(fp)]', beh(x).onRest, beh(x).offRest, 1); % identify sample during locomotion
    fp_imm = fp(idx_imm);

    %%
    a_on = align_full{x,y}([find(time == 0):find(time == 0.5)],:) - fp_mu; % SUBTRACT baseline%
    base = nanmean(align_full{x,y}([find(time == -4):find(time == -2)],:),1);
    a_on = a_on - base;
    a_on = a_on(:,all(~isnan(a_on))); % remove NaN trials
    [m, i] = max(a_on); i = i+find(time == 0)-1; % onset window: [0 1]s
    peak_onset{x} = m(:); % STORE
    
    prc_captured(x) = 0; counter = 0;
    while prc_captured(x) < 0.9 % Adjust threshold parameters until achieve 90%
        if counter >= 1; thres = thres - 0.2; end
        if thres < 1; break; end
        idx_onset = extractEventST([1:length(fp)]', beh(x).on, beh(x).on+25, 1); % Determine indices for onset window: [0 +0.5]
        fp_onset = fp(idx_onset); % Extract photometry signal during onset window: [0 +0.5]
        cross = find(fp_onset > thres); % Indices of samples where signal exceedes threshold
        [cross_cell, cross_pts] = consecutive_vec2cell(cross); % Extract grouped crossings
        peak_tmp = []; peak_tmp_idx = [];
        for c = 1:length(cross_cell) % Identify maximum of threshold crossing groups
            [peak_tmp(c),ii] = max(fp_onset(cross_cell{c}));
            peak_tmp_idx(c) = cross_cell{c}(ii);
        end 
        where = ismember(peak_tmp(:), m(:)+base); % Where is max onset
        peak_group = cross_cell(where); peak_pts = cross_pts(where,:); % Extract grouped crossings that match maximum onset group
        prc_captured(x) = length(find(where))/size(a_on,2); % Percent of onsets captured by threshold
        % peak_onset{x} = peak_tmp(where)'; % STORE only peaks that are above threshold
        counter = counter + 1; % Add to counter
    end
    thres_keep(x) = thres;
    %%
    % [idxPeakMax, idxPeakStart, peakMag, crossGroups, crossStart] = findIdxPeak(fpInputVec, thres, width)
    [idx_c_max, idx_c, peak_non{x}, ~, cross_pts] = findIdxPause(fp_imm, 'peak', thres, width); 

%%
%% ach to peaks
    idxpeak_imm = findIdxPause(fp(idx_imm), 'peak', thres, width);
    % idxpeak_mov = findIdxPause(fp(idx_onset), 'peak', thres, width);
    [ach2achpeak{x,1},sta_time] = getSTA(fp(idx_imm), idxpeak_imm/Fs, Fs, [-2 2]);
    % ach2achpeak{x,2} = getSTA(fp(idx_onset), idxpeak_mov/Fs, Fs, [-2 2]);
    ach2achpeak{x,2} = getSTA(fp, idx_onset(peak_tmp_idx)/Fs, Fs, [-2 2]);
    
    vel = beh(x).vel; acc = getAcc(vel);
    vel2achpeak{x,1} = getSTA(vel(idx_imm), idxpeak_imm/Fs, Fs, [-2 2]);
    vel2achpeak{x,2} = getSTA(vel, idx_onset(peak_tmp_idx)/Fs, Fs, [-2 2]);
    
    acc2achpeak{x,1} = getSTA(acc(idx_imm), idxpeak_imm/Fs, Fs, [-2 2]);
    acc2achpeak{x,2} = getSTA(acc, idx_onset(peak_tmp_idx)/Fs, Fs, [-2 2]);
    
%% 
waitbar(x/length(beh),h);
end
close(h);

% n = X mice
tmp = {}; for x = 1:length(beh); tmp{x} = strtok(beh(x).rec,'-'); end
uni = unique(tmp); nAn = length(uni);
ach2peak_an = cell(2,1); vel2peak_an = cell(2,1); acc2peak_an = cell(2,1);

for x = 1:nAn
    ii = find(strcmp(tmp,uni{x}));
    for z = 1:2
        ach2peak_an{z}(:,x) = nanmean([ach2achpeak{ii,z}],2);
        vel2peak_an{z}(:,x) = nanmean([vel2achpeak{ii,z}],2);
        acc2peak_an{z}(:,x) = nanmean([acc2achpeak{ii,z}],2);
    end
end

%
fig = figure; fig.Position(3) = 1375;

subplot(1,3,1); hold on; clr = {'k','g'};
plot([0 0],[-2 12],'--k');
for y = 1:length(ach2peak_an)
shadederrbar(sta_time, nanmean(ach2peak_an{y},2), SEM(ach2peak_an{y},2), clr{y}); hold on
end
xlabel('Latency to ACh peak (s)'); xlim([-1 1]); xticks([-1:0.5:1]);
ylabel('ACh3.0 (%dF/F)'); yticks([-4:4:14]);
axis('square')
title(sprintf('ACh to ACh peak (n = %d)',nAn));

subplot(1,3,2); hold on; 
plot([0 0],[-2 10],'--k');
for y = 1:length(vel2peak_an)
shadederrbar(sta_time, nanmean(vel2peak_an{y},2), SEM(vel2peak_an{y},2), clr{y}); hold on
end
xlabel('Latency to ACh peak (s)'); xlim([-1 1]); xticks([-1:0.5:1]);
ylabel('Velocity (cm/s)'); % yticks([-4:4:14]);
axis('square')
title(sprintf('Velocity to ACh peak'));

subplot(1,3,3); hold on;
plot([0 0],[-0.1 0.5],'--k');
for y = 1:length(acc2peak_an)
shadederrbar(sta_time, nanmean(acc2peak_an{y},2), SEM(acc2peak_an{y},2), clr{y}); hold on
end
xlabel('Latency to ACh peak (s)'); xlim([-1 1]); xticks([-1:0.5:1]);
ylabel('Acceleration (cm/s^2)');
axis('square')
title(sprintf('Acceleration to ACh peak'));
movegui(fig,'center');