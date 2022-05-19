%% summary
% Find instantaneous phase of each signal
% Plot distribution of acceleration w.r.t. instantaneous phase:
% UPDATED: 16-May-2022

%%
a = [1:length(modACh)]; rmv = [3 17 20 22 23 30 36:39 41:43 48 50]; a(rmv) = [];
beh = modACh(a); % Extract recordings with reward

%%
fp2ph_beh = cell(2,3);
fp2ph_prc_beh = cell(3,1);

h = waitbar(0, 'photometry to instantaneous phase');
for x = 1:length(beh)
    
    acc = getAcc(beh(x).vel); % Acceleration
    
    fp_mat = []; fp_phase = []; fp_deg = [];
    fp2ph = cell(2,2); for a = 1:2; for b = 1:2; fp2ph{a,b} = nan(36, 10000); end; end % initialize matrix
    fp2ph_state = cell(2,3); for a = 1:2; for b = 1:3; fp2ph_state{a,b} = nan(36, 10000); end; end % initialize matrix

    for y = 1
    fp_mat(:,y) = beh(x).FP{y} - nanmean(beh(x).FP{y}); % Photometry

    %% Index of behavioral states
    idx_cell = cell(3,1); 
    if isfield(beh,'reward')
        idx_rew = extractEventST([1:length(fp_mat)]', floor(beh(x).reward), floor(beh(x).reward)+100, 1); % identify sample during reward
        idx_cell{3} = idx_rew; 
    else; idx_rew = []; idx_cell = cell(2,1); end
    idx_mov = extractEventST([1:length(fp_mat)]', beh(x).on, beh(x).off, 1); % identify sample during locomotion
    idx_mov_nonRew = idx_mov(~ismember(idx_mov, idx_rew)); % exclude reward, include locomotion
    idx_imm = extractEventST([1:length(fp_mat)]', beh(x).onRest, beh(x).offRest, 1); % identify sample during locomotion
    idx_imm_nonRew = idx_imm(~ismember(idx_imm, idx_rew)); % exclude reward, include rest
    idx_cell{1} = idx_imm_nonRew; idx_cell{2} = idx_mov_nonRew; % index into cell array for ease of iteration

    %% Bandpass filter
    signal = fp_mat(:,y);
    Fpass = [0.5 4];
    Fs = 50; %sampling rate, has to be at least double of your high pass frequency
    Wn = [Fpass(1)/(Fs/2) Fpass(2)/(Fs/2)]; 
    [b,a] = butter(3,Wn);
    data_filt= filtfilt(b,a,signal); % signal is your lfp data, output is the filtered data

    %% Instantaneous phase
    H = hilbert(double(data_filt));
    data_phase = angle(H); % output is the instantaneous phase
    fp_phase(:,y) = data_phase;
    fp_deg(:,y) = rad2deg(data_phase);

    %Take the absolute of the filtered signal and calculate the standard deviation
    rmssig  = abs(data_filt);
    stdsig = std(rmssig);

    %Find the index of peaks that are more than NumStd standard deviations (NumStd could be 1.5 standard deviations)
    % finds the 0 degree phase indices that cross 1.5 standard deviations. Might have to remove the ; depending on whether your data is a row or column array
    NumStd = 1.5;
    peakIdx = find(data_filt>(NumStd*stdsig) & [0; diff(data_phase>0)]);         
        
    %% Align ACCELERATION to phase
    [~, edges] = histcounts(fp_deg(:,y), 'BinWidth', 10, 'BinLimits', [-180 180]); % genereate edges vector for phase, in degrees
    [phIdx] = discretize(fp_deg(:,y), edges); % returns the indices of the bins that the elements of X fall into

    for z = 1:length(unique(phIdx))
        ii = find(phIdx == z); % find index in signal that are this phase
        fp2ph{y,1}(z,[1:length(ii)]) = acc(ii); % extract ACCELERATION values matching this phase
        for c = 1:2
            ii_c = ii(ismember(ii, idx_cell{c})); % include index for this behavioral state
            fp2ph_state{y,c}(z,[1:length(ii_c)]) = acc(ii_c); % extract ACCELERATION values for this behavioral state
        end
        c = 3;
        acc_rand = acc(randperm(length(acc),length(acc)));
        fp2ph_state{y,c}(z,[1:length(ii_c)]) = acc_rand(ii_c); % extract ACCELERATION values for these random indices
    end
    
    [~,b] = find(isnan(fp2ph{y,1})); % find first column with NaN's
    fp2ph{y,1}(:,[b(1):size(fp2ph{y,1},2)]) = []; % remove end of matrix with incomplete columns
    fp2ph{y,2} = normalize(fp2ph{y,1},1,'range'); % normalize over each individual [0 180] oscillation
    % [~,b] = find(isnan(fp2ph{y,2})); % find first column with NaN's
    % fp2ph{y,2}(:,[b(1):size(fp2ph{y,2},2)]) = []; % remove end of matrix with incomplete columns
    
    %% save into output cell array
    for z = 1:length(idx_cell)
        fp2ph_beh{y,z}(:,x) = nanmean(fp2ph_state{y,z},2);
    end
    fp2ph_prc_beh{1}(:,x) = prctile(fp2ph_state{y,3},25,2);
    fp2ph_prc_beh{2}(:,x) = prctile(fp2ph_state{y,3},50,2);
    fp2ph_prc_beh{3}(:,x) = prctile(fp2ph_state{y,3},75,2);
    end
    %%
    waitbar(x/length(beh),h);
end
fprintf('Done! \n'); close(h);
% clearvars -except modAChDA beh edges fp2ph_beh fp2ph_norm_beh da2achph_norm_beh

% N = X mice: DA to ACh pause
tmp = {}; for x = 1:length(beh); tmp{x} = strtok(beh(x).rec,'-'); end
uni = unique(tmp); nAn = length(uni);
fp2ph_an = cell(2,3); % Initialize cell arrays
fp2ph_prc_an = cell(3,1);

for x = 1:nAn
    ii = find(strcmp(tmp,uni{x}));
    for z = 1:length(idx_cell)
        for y = 1
            fp2ph_an{y,z}(:,x) = nanmean([fp2ph_beh{y,z}(:,ii)],2);
            fp2ph_prc_an{z}(:,x) = nanmean(fp2ph_prc_beh{z}(:,ii),2);
        end
    end
end

%% PLOT
fig = figure; hold on; 
clr = {'r','g','b'};
mid = edges(2:end)-5;
y = 1; switch y; case 1; lbl = 'ACh'; case 2; lbl = 'DA'; end
for z = 1:2
    shadederrbar( mid, nanmean(fp2ph_an{y,z},2), SEM(fp2ph_an{y,z},2), clr{z});
end
a = nanmean(fp2ph_prc_an{2} - fp2ph_prc_an{1},2);
b = nanmean(fp2ph_prc_an{3} - fp2ph_prc_an{2},2);
shadederrbar(mid, nanmean(fp2ph_prc_an{2},2), nanmean(fp2ph_prc_an{1},2), 'k');

ylabel('Acceleration (cm/s^2)'); %ylim([-0.1 0.15]);
xlabel(sprintf('%s phase (deg)',lbl)); xlim([-180 180]); xticks([-180:45:180]);
title(sprintf('Acceleration @ %s phase',lbl)); axis('square');

%% QUANTIFICATION
pref = cell(2,1);

for a = 1:size(fp2ph_an,1); for b = 1:size(fp2ph_an,2)
    [~,c] = max(fp2ph_an{a,b}); % maximum, normalized to 1 so this is the preferred phase
    pref{a}(:,b) = mid(c)+180; % preferred phase at maximum
    end; end

fig = figure; fig.Position(3) = 1000; clr = {'r','g','b'};
for a = 1:2
    subplot(1,2,a,polaraxes); hold on
    for b = 1:2
    polarhistogram(deg2rad(pref{a}(:,b)),6,'FaceColor',clr{b},'FaceAlpha',0.3,'Normalization','count');
    end
end