%%
thres_pause = -3;
width_pause = 5;
thres_peak = 5; % CHANGE -- threshold for peaks/pauses to exceed, after baseline
width_peak = 5; % CHANGE -- require that peaks/pauses are at least X samples long (1 sample = 20ms)

h = waitbar(0, 'find immobility pauses/peaks');
for y = 1:length(s) % iterate over infusion
    %%
    beh = s(y).s;
    amp = cell(length(beh),2); 
    dur = cell(length(beh),2); freq = [];
    ach2ach = cell(length(beh),2);
    da2ach = cell(length(beh),2);
    params = [];

    for x = 1:length(beh) % iterate over animal
        if isempty(beh(x).Fs); continue; end
        fp = beh(x).FP{1}; Fs = 50; % extract full photometry signal
        idx_inf = [s(y).win(x,1).*(50*60) : s(y).win(x,2).*(50*60)]'; % infusion window
        rewWindow = 100; % how many samples after reward delivery is the reward window
        idx_inf_rew = extractEventST(idx_inf, floor(beh(x).reward), floor(beh(x).reward)+rewWindow, 1); % identify sample during reward
        idx_inf_rest = extractEventST(idx_inf, beh(x).onRest, beh(x).offRest, 1); % index: infusion + rest
        idx_inf_rest = idx_inf_rest(~ismember(idx_inf_rest, idx_inf_rew)); % exclude reward, include rest
        % fp = fp - s(y).base_mu_rest(x); % subtract mean of baseline recording to center around 0%
        fp = fp - nanmean(fp);
        
    %% pause identification
        switch strtok(beh(x).rec,'_')
            case 'AK189'; thres_peak = 2.75; case 'AK190'; thres_peak = 3; case 'AK193'; thres_peak = 4.5; case 'AK197'; thres_peak = 4; 
        end
        switch strtok(beh(x).rec,'_')
            case 'AK189'; thres_pause = -2.25; case 'AK190'; thres_pause = -2.5; case 'AK193'; thres_pause = -2; case 'AK197'; thres_pause = -2.25; 
        end
    
        [pause_maxIdx, pause_good, pause_valmag, ~, pause_crossstart] = findIdxPause(fp(idx_inf_rest),'pause',thres_pause,width_pause);
        [peak_maxIdx, peak_good, peak_valmag, ~, peak_crossstart] = findIdxPause(fp(idx_inf_rew),'peak',thres_peak,width_peak);
        if y == 3
            if strtok(beh(x).rec,'_') == 'AK190'
                [pause_maxIdx, pause_good, pause_valmag, ~, pause_crossstart] = findIdxPause(fp(idx_inf_rest),'pause',-3.5,width_pause);
            end
        end
        if y == 4
            % [peak_maxIdx, peak_good, peak_valmag, ~, peak_crossstart] = findIdxPause(fp(idx_inf_rew),'peak',5,width_peak);
%             if x == 1
%                 [pause_maxIdx, pause_good, pause_valmag, ~, pause_crossstart] = findIdxPause(fp(idx_inf_rest),'pause',-2.25,width_pause);
%             else
%                 [pause_maxIdx, pause_good, pause_valmag, ~, pause_crossstart] = findIdxPause(fp(idx_inf_rest),'pause',-3.5,width_pause);
%             end
        end
        params(x,:) = [thres_pause, thres_peak];
        
    %% ACh to ACh pause
        [ach2ach{x,1},sta_time] = getSTA(fp(idx_inf_rest), pause_maxIdx/Fs, Fs, [-6 2]);
        ach2ach{x,2} = getSTA(fp(idx_inf_rew), peak_maxIdx/Fs, Fs, [-6 2]);
        
        fp_da = beh(x).FP{2};
        fp_da = fp_da - nanmean(fp_da);
        da2ach{x,1} = getSTA(fp_da(idx_inf_rest), pause_maxIdx/Fs, Fs, [-6 2]);
        da2ach{x,2} = getSTA(fp_da(idx_inf_rew), peak_maxIdx/Fs, Fs, [-6 2]);
        
    %%
        amp{x,1} = pause_valmag(pause_good); % Pause magnitude
        amp{x,2} = peak_valmag(peak_good); % Peak magnitude
        dur{x,1} = diff(pause_crossstart(pause_good,:),1,2); % Duration of pause that exceeds threshold
        dur{x,2} = diff(peak_crossstart(peak_good,:),1,2); % Duration of pause that exceeds threshold
        freq(x,1) = (1/mean(diff(pause_maxIdx)))*Fs; % Frequency of maximum deflections
        freq(x,2) = (1/mean(diff(peak_maxIdx)))*Fs; % Frequency of maximum deflections
        
    waitbar(x/length(beh),h);
%%
    end
    s(y).thres = params; s(y).width = [width_pause, width_peak];
    s(y).amp = amp; 
    s(y).dur = dur; 
    s(y).freq = freq;
    s(y).ach2ach = ach2ach;
    s(y).da2ach = da2ach;
    
end
close(h);