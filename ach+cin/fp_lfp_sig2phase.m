%% summary
% Find instantaneous phase of each signal
% Plot distribution of %dF/F photometry signal w.r.t. instantaneous phase:
%   for ACh photometry signal, normalized, to the phase of LFP

%%
load('C:\Users\Anya\Desktop\IV_LOCAL\data_comb\fullbeh_ACh_v2.mat')
beh = behACh([1 3 8 11]);
lfp_path = 'C:\Users\Anya\Desktop\IV_LOCAL\ACh';
[fName,fPath] = uigetfile([lfp_path],'MultiSelect','On');
for x = 1:length(beh)
    load(fullfile(fPath,fName{x}));
    if ~all(lfp.rec == beh(x).rec); error('fp and lfp do not match recording'); end
        
    lfpds = downsampleTLab(double(lfp.data), lfp.samplingRate/beh(x).Fs, 2); % downsample LFP signal to match sampling rate of photometry signal
    data_mat = beh(x).FP{1} - nanmean(beh(x).FP{1}); 
    data_mat(:,2) = lfpds(1:length(data_mat));
    beh(x).FP{2} = data_mat(:,2); beh(x).FPnames{2} = ['LFP-ch',num2str(lfp.channels)];
end

%%
ach_lfp = cell(2,2);
ach_norm_lfp = ach_lfp;

h = waitbar(0, 'photometry to instantaneous phase');
for x = 1:length(beh)
    % load(fullfile(fPath,fName{x}));
    % if ~all(lfp.rec == beh(x).rec); error('fp and lfp do not match recording'); end 
    % lfpds = downsampleTLab(double(lfp.data), lfp.samplingRate/beh(x).Fs, 2); % downsample LFP signal to match sampling rate of photometry signal
    data_mat = beh(x).FP{1} - nanmean(beh(x).FP{1}); % PHOTOMETRY
    data_mat(:,2) = beh(x).FP{2}; % LFP
    
    %% Index of behavioral states
    if isfield(beh,'reward')
        idx_rew = extractEventST([1:length(data_mat)]', floor(beh(x).reward), floor(beh(x).reward)+100, 1); % identify sample during reward
    else; idx_rew = []; end
    idx_mov = extractEventST([1:length(data_mat)]', beh(x).on, beh(x).off, 1); % identify sample during locomotion
    idx_mov_nonRew = idx_mov(~ismember(idx_mov, idx_rew)); % exclude reward, include locomotion
    idx_imm = extractEventST([1:length(data_mat)]', beh(x).onRest, beh(x).offRest, 1); % identify sample during locomotion
    idx_imm_nonRew = idx_imm(~ismember(idx_imm, idx_rew)); % exclude reward, include rest
    idx_cell = cell(3,1); idx_cell{1} = idx_imm_nonRew; idx_cell{2} = idx_mov_nonRew; idx_cell{3} = idx_rew; % index into cell array for ease of iteration

    %%
    data_peakIdx = cell(2,1);
    data_filt = []; data_deg = [];
    
    for y = 1:2
        %% Bandpass filter
        signal = data_mat(:,y);
        Fpass = [0.5 4];
        Fs = 50; %sampling rate, has to be at least double of your high pass frequency
        Wn = [Fpass(1)/(Fs/2) Fpass(2)/(Fs/2)]; 
        [b,a] = butter(3,Wn);
        data_filt(:,y) = filtfilt(b,a,signal); % signal is your lfp data, output is the filtered data

        %% Instantaneous phase
        H = hilbert(double(data_filt(:,y)));
        data_phase = angle(H); % output is the instantaneous phase
        data_deg(:,y) = rad2deg(data_phase);

        %Take the absolute of the filtered signal and calculate the standard deviation
        rmssig  = abs(data_filt(:,y));
        stdsig = std(rmssig);

        %Find the index of peaks that are more than NumStd standard deviations (NumStd could be 1.5 standard deviations)
        % finds the 0 degree phase indices that cross 1.5 standard deviations. Might have to remove the ; depending on whether your data is a row or column array
        NumStd = 1.5;
        peakIdx = find(data_filt(:,y)>(NumStd*stdsig) & [0; diff(data_phase>0)]);         
        data_peakIdx{y} = peakIdx;
    end
    
    %% Align one signal to phase of other
    for y = 1:2
        switch y; case 1; a = 1; b = 2; case 2; a = 2; b = 1; end
        [~, edges] = histcounts(data_deg(:,a), 'BinWidth', 10, 'BinLimits', [-180 180]); % genereate edges vector for phase, in degrees
        [phIdx] = discretize(data_deg(:,a), edges); % returns the indices of the bins that the elements of X fall into
        
        tmp_data = cell(1,3); for aa = 1:3; tmp_data{aa} = nan(36, 10000); end % initialize matrix
        tmp_norm = tmp_data;
        for z = 1:length(unique(phIdx))
            ii = find(phIdx == z); % find index in signal that are this phase
            for c = 1:2
                ii_c = ii(ismember(ii, idx_cell{c})); % include index for this behavioral state
                tmp_data{c}(z,[1:length(ii_c)]) = data_mat(ii_c,b); % extract photometry values for this behavioral state
            end
        end
        for z = 1:2
            [~,ii] = find(isnan(tmp_data{z})); % find first column with NaN's
            tmp_data{z}(:,[ii(1):size(tmp_data{z},2)]) = []; % remove end of matrix with incomplete columns
            ach_lfp{y,z}(:,x) = nanmean(tmp_data{z},2); % save into output cell array
            
            tmp_norm{z} = normalize(tmp_data{z},1,'range');
%             [~,ii] = find(isnan(tmp_norm{z})); % find first column with NaN's
%             tmp_norm{z}(:,[ii(1):size(tmp_norm{z},2)]) = []; % remove end of matrix with incomplete columns
            ach_norm_lfp{y,z}(:,x) = normalize(nanmean(tmp_norm{z},2),1,'range'); % save into output cell array
        end
        
    end
    %%
    waitbar(x/length(beh),h);
end
mid = edges(2:end)'-5; % midpoints of edge bins, for plotting
fprintf('Done! \n'); close(h);

%% PLOT IMMOBILITY normalized signal to phase
fig = figure; fig.Position(3) = 1375; movegui(gcf, 'center')
sm = 5;

for y = 1:2
    % switch y; case 1; lbl = {'ACh','LFP'}; case 2; lbl = {'LFP','ACh'}; end
    switch y; case 1; lbl = {'DA','LFP'}; case 2; lbl = {'LFP','DA'}; end

    subplot(1,3,y); hold on
    z = 1; clr = 'm'; % clr = [0.05 0.75 0.45];
    a = ach_norm_lfp{y,z};
    b = [a(19:36,:);a(1:18,:)];
    
    plot( mid+180, movmean(b,sm,1), 'Color', [0 0 0 0.5]);
    plot( mid+540, movmean(b,sm,1), 'Color', [0 0 0 0.5]);
    shadederrbar( mid+180, movmean(nanmean(b,2),sm), movmean(SEM(b,2),sm), clr);
    shadederrbar( mid+540, movmean(nanmean(b,2),sm), movmean(SEM(b,2),sm), clr);

    xlabel(sprintf('%s Phase',lbl{1})); xlim([0 720]); xticks([0:180:720]);
    ylabel(sprintf('%s Signal (norm.)',lbl{2})); ylim([0 1]);
    title(sprintf('IMM %s norm to %s phase',lbl{2},lbl{1})); axis('square');
end

a = ach_norm_lfp{2,1}; % extract for ACh to LFP phase
b = [a(19:36,:);a(1:18,:)];
[c,~] = find(b == 1); % find maximum
pref = mid(c)+180;
subplot(1,3,3); % Phase preferences for all units
rose(deg2rad(pref));
title('IMM fp preferred phase of LFP')

movegui(gcf,'center');

%% PLOT by ANIMAL -- norm and raw
fig = figure; fig.Position([3 4]) = [1000 850]; movegui(gcf, 'center')
clr = {'r','g'};
sm = 5;
for y = 1:2
    switch y; case 1; lbl = {'ACh','LFP'}; case 2; lbl = {'LFP','ACh'}; end
    
    a = y; sp(y) = subplot(2,2,a); hold on
    for z = 1:2
        shadederrbar( mid, movmean(nanmean(ach_lfp{y,z},2),sm), movmean(SEM(ach_lfp{y,z},2),sm), clr{z});
    end
    xlabel(sprintf('%s Phase',lbl{1})); xlim([-180 180]); xticks([-180:90:180]);
    ylabel(sprintf('%s Signal',lbl{2}));
    % legend({'imm','','mov'});
    title(sprintf('%s signal to %s phase (n = %d)',lbl{2},lbl{1},size(ach_lfp{1,1},2))); axis('square');
    
    a = y+2; subplot(2,2,a); hold on
    for z = 1:2
        shadederrbar( mid, movmean(nanmean(ach_norm_lfp{y,z},2),sm), movmean(SEM(ach_norm_lfp{y,z},2),sm), clr{z});
        %plot( mid, movmean(nanmean(ach_norm_lfp{y,z},2),sm), clr{z});
    end
    xlabel(sprintf('%s Phase',lbl{1})); xlim([-180 180]); xticks([-180:90:180]);
    ylabel(sprintf('%s Signal (norm.)',lbl{2})); ylim([0 1]);
    % legend({'imm','','mov'});
    title(sprintf('%s norm to %s phase',lbl{2},lbl{1})); axis('square');
end
