%% summary
% For one signal, filter between 0.5-4 Hz and detect peaks
% If you have the instantaneous phase of each signal
%   Calculate the average phase for signal 1 for about 500 ms on either side of the peak events(like a PETH). 
%   Then calculate the average phase for signal 2 using the detected peaks for signal 1.  


%% 
ph_da = cell(length(beh),3); ph_ach = cell(length(beh),3); 
for x = 1:length(beh)
    fp_mat = []; fp_phase = []; fp_peak = cell(2,3); fp_avgph = cell(2,1);
    for y = 1:2
        fp_mat(:,y) = beh(x).FP{y} - nanmean(beh(x).FP{y}); 

        %% Index of behavioral states
        if isfield(beh,'reward')
            idx_rew = extractEventST([1:length(fp_mat)]', floor(beh(x).reward)-100, floor(beh(x).reward)+100, 1); % identify sample during reward
        else; idx_rew = []; end
        idx_mov = extractEventST([1:length(fp_mat)]', beh(x).on, beh(x).off, 1); % identify sample during locomotion
        idx_mov_nonRew = idx_mov(~ismember(idx_mov, idx_rew)); % exclude reward, include locomotion
        idx_imm = extractEventST([1:length(fp_mat)]', beh(x).onRest, beh(x).offRest, 1); % identify sample during locomotion
        idx_imm_nonRew = idx_imm(~ismember(idx_imm, idx_rew)); % exclude reward, include rest
        c = cell(3,1); c{1} = idx_imm_nonRew; c{2} = idx_mov_nonRew; c{3} = idx_rew; % index into cell array for ease of iteration

        %% Bandpass filter
        signal = fp_mat(:,y);
        Fpass = [0.5 4];
        Fs = 50; %sampling rate, has to be at least double of your high pass frequency
        Wn = [Fpass(1)/(Fs/2) Fpass(2)/(Fs/2)]; 
        [b,a] = butter(3,Wn);
        data_filt= filtfilt(b,a,signal); % signal is your lfp data, output is the filtered data

        %% Instantaneous phase
        H = hilbert(double(data_filt));
        data_phase  = angle(H); % output is the instantaneous phase
        fp_phase(:,y) = data_phase;

        %Take the absolute of the filtered signal and calculate the standard deviation
        rmssig  = abs(data_filt);
        stdsig = std(rmssig);

        %Find the index of peaks that are more than NumStd standard deviations (NumStd could be 1.5 standard deviations)
        % finds the 0 degree phase indices that cross 1.5 standard deviations. Might have to remove the ; depending on whether your data is a row or column array
        NumStd = 1.5;
        peakIdx = find(data_filt>(NumStd*stdsig) & [0; diff(data_phase>0)]);
        for z = 1:3
            peakIdx_sub = peakIdx(ismember(peakIdx,c{z}));
            fp_peak{y,z} = peakIdx_sub;
        end
    end
    y = 1; % instantaneous phase at ACh(y=1) peaks
    for z = 1:3 
        ph_ach{x,z} = rad2deg(fp_phase(fp_peak{y,z},1));
        ph_da{x,z} = rad2deg(fp_phase(fp_peak{y,z},2)); end
end

phase_cell = cell(2,3);
for x = 1:size(ph_da,1)
    for z = 1:3
        [phase_cell{1,z}(:,x), edges] = histcounts(ph_ach{x,z},'BinWidth',10,'BinLimits',[-180 180],'Normalization','probability');
        [phase_cell{2,z}(:,x), edges] = histcounts(ph_da{x,z},'BinWidth',10,'BinLimits',[-180 180],'Normalization','probability');
    end 
end

% N = X mice
phase_an = cell(2,3);
tmp = {}; for x = 1:length(beh); tmp{x} = strtok(beh(x).rec,'-'); end
uni = unique(tmp); nAn = length(uni);

for x = 1:nAn
    idx = strcmp(tmp,uni{x});
    for a = 1:2; for b = 1:3
        phase_an{a,b}(:,x) = nanmean(phase_cell{a,b}(:,idx),2);
        end; end
end

fprintf('Done.\n');
%%
y = 2; %ACh inst phase @ DA peak

fig = figure; fig.Position(3) = 1400;
mid = edges(2:end)-5;
subplot(1,3,1);
a = phase_an{y,1};
c = a'; % m = max(a); [~,b] = sort(m); c = a(:,b)';
z = 1; h(z) = heatmap(c); h(z).YLabel = 'recording';
h(z).Title = 'IMM only: ACh inst phase @ DA peak';
h(z).XLabel = 'phase';
h(z).Colormap = jet; h(z).GridVisible = 'off'; % h(z).ColorLimits = [0 1];

subplot(1,3,2); hold on
shadederrbar(mid, nanmean(phase_an{y,1},2), SEM(phase_an{y,1},2), 'g');
yyaxis right; shadederrbar(mid, nanmean(phase_an{2,1},2), SEM(phase_an{2,1},2), 'm');
legend({'ACh','','DA'},'Location','northwest');
xlabel('phase'); xlim([-180 180]); xticks([-180:90:180]); 
title('ACh inst phase @ DA peak'); axis('square');

sm = 5;
subplot(1,3,3); hold on
shadederrbar(mid, movmean(nanmean(phase_an{y,1},2),sm), movmean(SEM(phase_an{y,1},2),sm), 'r');
shadederrbar(mid, movmean(nanmean(phase_an{y,2},2),sm), movmean(SEM(phase_an{y,2},2),sm), 'g');
shadederrbar(mid, movmean(nanmean(phase_an{y,3},2),sm), movmean(SEM(phase_an{y,3},2),sm), 'b');
legend({'imm','','mov','','rew'},'Location','northwest');
xlabel('phase'); xlim([-180 180]); xticks([-180:90:180]); 
title('ACh inst phase @ DA peak'); axis('square');

%% TESTING
x = 1;

fp_mat = []; fp_phase = []; fp_peak = cell(2,1); fp_avgph = cell(2,1);
for y = 1:2
    fp_mat(:,y) = beh(x).FP{y};
    fp_mat(:,y) = fp_mat(:,y) - nanmean(beh(x).FP{y}); 
    
    %% Index of behavioral states
    if isfield(beh,'reward')
        idx_rew = extractEventST([1:length(fp_mat)]', floor(beh(x).reward)-100, floor(beh(x).reward)+100, 1); % identify sample during reward
    else; idx_rew = []; end
    idx_mov = extractEventST([1:length(fp_mat)]', beh(x).on, beh(x).off, 1); % identify sample during locomotion
    idx_mov_nonRew = idx_mov(~ismember(idx_mov, idx_rew)); % exclude reward, include locomotion
    idx_imm = extractEventST([1:length(fp_mat)]', beh(x).onRest, beh(x).offRest, 1); % identify sample during locomotion
    idx_imm_nonRew = idx_imm(~ismember(idx_imm, idx_rew)); % exclude reward, include rest
    c = cell(3,1); c{1} = idx_imm_nonRew; c{2} = idx_mov_nonRew; c{3} = idx_rew; % index into cell array for ease of iteration
    
    %% Bandpass filter
    signal = fp_mat(:,y);
    Fpass = [0.5 4];
    Fs = 50; %sampling rate, has to be at least double of your high pass frequency
    Wn = [Fpass(1)/(Fs/2) Fpass(2)/(Fs/2)]; 
    [b,a] = butter(3,Wn);
    data_filt= filtfilt(b,a,signal); % signal is your lfp data, output is the filtered data

    %% Instantaneous phase
    H = hilbert(double(data_filt));
    data_phase  = angle(H); % output is the instantaneous phase
    fp_phase(:,y) = data_phase;
    
    %Take the absolute of the filtered signal and calculate the standard deviation
    rmssig  = abs(data_filt);
    stdsig = std(rmssig);

    %Find the index of peaks that are more than NumStd standard deviations (NumStd could be 1.5 standard deviations)
    % finds the 0 degree phase indices that cross 1.5 standard deviations. Might have to remove the ; depending on whether your data is a row or column array
    NumStd = 1.5;
    peakIdx = find(data_filt>(NumStd*stdsig) & [0; diff(data_phase>0)]);         
    fp_peak{y} = peakIdx;
    
    %% average phase for signal on either side of the peak events
    % [fp_avgph{y}, t] = getSTA(rad2deg(fp_phase(:,y)), fp_peak{1}./Fs, Fs, [-2 2]);
end

%%
figure; hold on
histogram(rad2deg(fp_phase(fp_peak{1},1)),'BinWidth',10,'Normalization','probability','FaceColor','g','FaceAlpha',0.2);
yyaxis right; histogram(rad2deg(fp_phase(fp_peak{1},2)),'BinWidth',10,'Normalization','probability','FaceColor','m','FaceAlpha',0.2);
legend({'ACh @ ACh peak','DA @ ACh peak'});
xlabel('phase'); xlim([-180 180]); xticks([-180:90:180]);
title('inst phase at ACh peak')