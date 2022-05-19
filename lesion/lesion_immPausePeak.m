%% summary
%%
for y = 1:2
    %%
    beh = s(y).s; % extract beh sub-structure within cannula structure
    amp = cell(length(beh),2); 
    dur = cell(length(beh),2); 
    freq = []; params = [];
    ach2ach = cell(length(beh),2);
    da2ach = cell(length(beh),2);

    for x = 1:length(beh) % iterate over animal
        if isempty(beh(x).Fs); continue; end
        fp_mat = [beh(x).FP{1}, beh(x).FP{2}]; Fs = 50;
        fp_mat = fp_mat - nanmean(fp_mat);
                
        rewWindow = Fs; % how many samples after reward delivery is the reward window
        idx_rew = extractEventST([1:length(fp_mat)]', floor(beh(x).reward), floor(beh(x).reward)+rewWindow, 1); % identify sample during reward
        idx_rest = extractEventST([1:length(fp_mat)]', beh(x).onRest, beh(x).offRest, 1); % index: infusion + rest
        idx_rest = idx_rest(~ismember(idx_rest, idx_rew)); % exclude reward, include rest

        %% Bandpass filter
        Fpass = [0.5 4];
        Fs = 50; %sampling rate, has to be at least double of your high pass frequency
        Wn = [Fpass(1)/(Fs/2) Fpass(2)/(Fs/2)]; 
        [b,a] = butter(3,Wn);
        data_filt = filtfilt(b, a, fp_mat(:,1)); % signal is your photometry data, output is the filtered data

        %% Instantaneous phase
        data_filt = data_filt(idx_rest); % IMMOBILITY ONLY
        H = hilbert(double(data_filt));
        data_phase  = angle(H); % output is the instantaneous phase
        fp_phase = data_phase;
        fp_deg = rad2deg(data_phase);

        %Take the absolute of the filtered signal and calculate the standard deviation
        rmssig  = abs(data_filt);
        stdsig = std(rmssig);        

        %Find the index of peaks that are more than NumStd standard deviations (NumStd could be 1.5 standard deviations)
        % finds the 0 degree phase indices that cross 1.5 standard deviations. Might have to remove the ; depending on whether your data is a row or column array
        NumStd = 2;
        if y == 1; params(x) = NumStd*stdsig;
        else
            ii = find(strcmp(uni,strtok(beh(x).rec,'-')));
            params(x) = s(1).params_an(ii); 
        end
            idxPeak = find(data_filt>(params(x)) & [0; diff(data_phase>0)]); 
            idxPause = find(-data_filt>(params(x)) & [0; diff(-data_phase>0)]); 
%         figure; hold on
%         plot(data_filt, 'k')
%         stem(idxPeak, 10*ones(length(idxPeak),1), 'g'); 
%         stem(idxPause, 10*ones(length(idxPause),1), 'r');

        %% Peak characterization
        tmp_amp_peak = data_filt(idxPeak); 
        tmp_amp_pause = data_filt(idxPause);
        tmp_dur_peak = []; tmp_dur_pause = []; win = 1*Fs;
        for z = 1:length(idxPeak)
            halfMax = 0.5*data_filt(idxPeak(z)); % amplitude at half-max
            if idxPeak(z)-win <= 0
                a = [data_filt(1 : idxPeak(z))]; a = flipud(a); % if peak is close to start of recording
            else
                a = [data_filt(idxPeak(z)-win : idxPeak(z))]; a = flipud(a); % segment preceding idx of maximum deflection
            end
            b = find(a < 0, 1, 'first') - 1; 
            if ~isempty(b); a = a(1:b); end % retain only values above 0
            c = []; [~,c(1)] = min(abs(a-halfMax)); c(1) = c(1)-1; % idx at half max
            if idxPeak(z)+win > length(data_filt)
                a = [data_filt(idxPeak(z) : end)]; 
            else
                a = [data_filt(idxPeak(z) : idxPeak(z)+win)]; % segment following idx of maximum deflection
            end
                
            b = find(a < 0, 1, 'first') - 1;
            if ~isempty(b); a = a(1:b); end % retain only values above 0
            [~,c(2)] = min(abs(a-halfMax)); % idx at half max
            
            tmp_dur_peak(z) = sum(c);
        end % width at half max for peaks
        for z = 1:length(idxPause)
            halfMax = 0.5*data_filt(idxPause(z)); % amplitude at half-max
            if idxPause(z)-win <= 0
                a = [data_filt(1 : idxPause(z))]; a = flipud(a); % if peak is close to start of recording
            else
                a = [data_filt(idxPause(z)-win : idxPause(z))]; a = flipud(a); % segment preceding idx of maximum deflection
            end
            b = find(a > 0, 1, 'first') - 1;
            if ~isempty(b); a = a(1:b); end % retain only values above 0
            c = []; [~,c(1)] = min(abs(a-halfMax)); c(1) = c(1)-1; % idx at half max
            if idxPause(z)+win > length(data_filt)
                a = [data_filt(idxPause(z) : end)]; 
            else
                a = [data_filt(idxPause(z) : idxPause(z)+win)]; % segment following idx of maximum deflection
            end
            b = find(a > 0, 1, 'first') - 1;
            if ~isempty(b); a = a(1:b); end % retain only values above 0
            [~,c(2)] = min(abs(a-halfMax)); % idx at half max
            
            tmp_dur_pause(z) = sum(c);
        end % width at half max for pauses

        amp{x,1} = tmp_amp_pause; amp{x,2} = tmp_amp_peak; % Amplitude of maximum deflection
        dur{x,1} = tmp_dur_pause; dur{x,2} = tmp_dur_peak; % Duration at half max
        freq(x,1) = (1./nanmean(diff(idxPause)))*Fs; % Frequency of maximum deflections
        freq(x,2) = (1./nanmean(diff(idxPeak)))*Fs; % Frequency of maximum deflections
        
        %% Photometry to peaks/pauses
%         [sta, t] = getSTA(fp_mat(idx_rest,1), idxPause/Fs, Fs, [-6 2]);
%         ach2ach{x,1} = sta;
%         sta = getSTA(fp_mat(idx_rest,1), idxPeak/Fs, Fs, [-6 2]);
%         ach2ach{x,2} = sta;
%         sta = getSTA(fp_mat(idx_rest,2), idxPause/Fs, Fs, [-6 2]);
%         da2ach{x,1} = sta;
%         sta = getSTA(fp_mat(idx_rest,2), idxPeak/Fs, Fs, [-6 2]);
%         da2ach{x,2} = sta;
    end
    s(y).params = params;
    s(y).amp = amp; 
    s(y).dur = dur; 
    s(y).freq = freq;
%     s(y).ach2ach = ach2ach;
%     s(y).da2ach = da2ach;
    
    %% n = X mice
    tmp = {}; for x = 1:length(beh); tmp{x} = strtok(beh(x).rec,'-'); end
    uni = unique(tmp); nAn = length(uni);
    s(y).amp_an = cell(nAn,2); s(y).dur_an = cell(nAn,2); s(y).freq_an = [];
%     s(y).ach2ach_an = cell(nAn,2); s(y).da2ach_an = cell(nAn,2); 
    s(y).params_an = [];
    for x = 1:nAn
        idx = find(strcmp(tmp,uni{x}));
        s(y).freq_an(x,:) = nanmean(s(y).freq(idx,:));
        s(y).params_an(x) = max(s(y).params(idx));
        for a = 1:length(idx)
        for z = 1:2
        s(y).amp_an{x,z} = [s(y).amp_an{x,z}(:); s(y).amp{idx(a),z}(:)];
        s(y).dur_an{x,z} = [s(y).dur_an{x,z}(:); s(y).dur{idx(a),z}(:)];
%         s(y).ach2ach_an{x,z} = [s(y).ach2ach_an{x,z}, s(y).ach2ach{idx(a),z}];
%         s(y).da2ach_an{x,z} = [s(y).da2ach_an{x,z}, s(y).da2ach{idx(a),z}];
        end
        end
    end
end
fprintf('lesion imm pause peak - done! \n');

%% PLOT
fig = figure; fig.Position([3 4]) = [1375 800];

for z = 1:2 % PAUSE(1) or PEAK(2)
switch z; case 1; lbl = 'Pause'; case 2; lbl = 'Peak'; end
nAn = 6;
%IMM - frequency
freq_mat = [];
for y = 1:length(s); freq_mat(:,y) = s(y).freq_an(:,z); end
freq_mat(isnan(freq_mat)) = 0;
a = freq_mat; %./nanmean(freq_mat(:,1));

subplot(2,3,1+(3*(z-1))); hold on
errorbar(abs(nanmean(a)'),abs(SEM(a,1)'),'.k','MarkerSize',20);
plot([1.15; 1.85].*ones(2,nAn),abs(a'),'.k','MarkerSize',20);
plot([1.15; 1.85].*ones(2,nAn),abs(a(:,[1 2])'),':r','MarkerSize',20);
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'pre','lesion'});
ylabel(sprintf('%s Frequency (Hz)',lbl)); ylim([0 2]); yticks([0:0.5:3]);
p = signrank(a(:,1),a(:,2));
title(sprintf('IMM: %s Freq (p = %1.3f)',lbl,p)); axis('square');

%IMM - duration
dur_mat = [];
for y = 1:length(s); for x = 1:nAn; dur_mat(x,y) = nanmean(s(y).dur_an{x,z}); end; end
dur_mat(isnan(dur_mat)) = 0;
dur_mat = dur_mat.*(1000/Fs); % Adjust from samples to ms
r = [1 2]; a = dur_mat(:,r);%./nanmean(dur_mat(:,1));

subplot(2,3,2+(3*(z-1))); hold on
errorbar(abs(nanmean(a)'),abs(SEM(a,1)'),'.k','MarkerSize',20);
plot([1.15; 1.85].*ones(2,nAn),abs(a'),'.k','MarkerSize',20);
plot([1.15; 1.85].*ones(2,nAn),abs(a(:,[1 2])'),':r','MarkerSize',20);
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'pre','lesion'});
ylabel(sprintf('%s Duration',lbl));  ylim([0 400]); yticks([0:100:500]);
p = signrank(a(:,1),a(:,2));
title(sprintf('IMM: %s Dur (p = %1.3f)',lbl,p)); axis('square');

%IMM - amplitude
amp_mat = [];
for y = 1:length(s); for x = 1:nAn; amp_mat(x,y) = nanmean(s(y).amp_an{x,z}); end; end
amp_mat(isnan(amp_mat)) = 0;
switch z; case 1; amp_mat = -1*amp_mat; end
r = [1 2]; a = amp_mat(:,r);%./nanmean(amp_mat(:,1));

subplot(2,3,3+(3*(z-1))); hold on
errorbar(abs(nanmean(a)'),abs(SEM(a,1)'),'.k','MarkerSize',20);
plot([1.15; 1.85].*ones(2,nAn),abs(a'),'.k','MarkerSize',20);
plot([1.15; 1.85].*ones(2,nAn),abs(a(:,[1 2])'),':r','MarkerSize',20);
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'pre','lesion'});
ylabel(sprintf('%s Amplitude',lbl)); ylim([0 5]); yticks([0:5]);
p = signrank(amp_mat(:,1),amp_mat(:,2));
title(sprintf('IMM: %s Amp (p = %1.3f)',lbl,p)); axis('square');
end
movegui(gcf,'center');