%% summary
% Find instantaneous phase of each signal
% Plot distribution of %dF/F photometry signal w.r.t. instantaneous phase:
%   (A) for each signal to its own phase, comparing behavioral states
%   (B) for each signal, normalized, to its own phase to demonstrate
%   relationship with focus on temporal aspect
%   (C) for DA photometry signal, normalized, to the phase of ACh

%%
x = 30;

fp_mat = []; fp_phase = []; fp_deg = [];
fp2ph = cell(2,2); for a = 1:2; for b = 1:2; fp2ph{a,b} = nan(36, 3000); end; end % initialize matrix
fp2ph_state = cell(2,3); for a = 1:2; for b = 1:3; fp2ph_state{a,b} = nan(36, 3000); end; end % initialize matrix
fp2ph_state_norm = cell(2,3); for a = 1:2; for b = 1:3; fp2ph_state_norm{a,b} = nan(36, 3000); end; end % initialize matrix

for y = 1:2
    fp_mat(:,y) = beh(x).FP{y} - nanmean(beh(x).FP{y}); 
    
    %% Index of behavioral states
    if isfield(beh,'reward')
        idx_rew = extractEventST([1:length(fp_mat)]', floor(beh(x).reward), floor(beh(x).reward)+100, 1); % identify sample during reward
    else; idx_rew = []; end
    idx_mov = extractEventST([1:length(fp_mat)]', beh(x).on, beh(x).off, 1); % identify sample during locomotion
    idx_mov_nonRew = idx_mov(~ismember(idx_mov, idx_rew)); % exclude reward, include locomotion
    idx_imm = extractEventST([1:length(fp_mat)]', beh(x).onRest, beh(x).offRest, 1); % identify sample during locomotion
    idx_imm_nonRew = idx_imm(~ismember(idx_imm, idx_rew)); % exclude reward, include rest
    idx_cell = cell(3,1); idx_cell{1} = idx_imm_nonRew; idx_cell{2} = idx_mov_nonRew; idx_cell{3} = idx_rew; % index into cell array for ease of iteration
    
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
    fp_deg(:,y) = rad2deg(data_phase);
    
    %Take the absolute of the filtered signal and calculate the standard deviation
    rmssig  = abs(data_filt);
    stdsig = std(rmssig);

    %Find the index of peaks that are more than NumStd standard deviations (NumStd could be 1.5 standard deviations)
    % finds the 0 degree phase indices that cross 1.5 standard deviations. Might have to remove the ; depending on whether your data is a row or column array
    NumStd = 1.5;
    peakIdx = find(data_filt>(NumStd*stdsig) & [0; diff(data_phase>0)]);         

    %% Align photometry to phase
    [~, edges] = histcounts(fp_deg(:,y), 'BinWidth', 10, 'BinLimits', [-180 180]); % genereate edges vector for phase, in degrees
    [phIdx] = discretize(fp_deg(:,y), edges); % returns the indices of the bins that the elements of X fall into

    for z = 1:length(unique(phIdx))
        ii = find(phIdx == z); % find index in signal that are this phase
        fp2ph{y,1}(z,[1:length(ii)]) = fp_mat(ii,y); % extract photometry %dF/F values matching this phase
        for c = 1:3
            ii_c = ii(ismember(ii, idx_cell{c})); % include index for this behavioral state
            fp2ph_state{y,c}(z,[1:length(ii_c)]) = fp_mat(ii_c,y); % extract photometry values for this behavioral state
        end
    end
    fp2ph{y,2} = normalize(fp2ph{y,1},1,'range'); % normalize over each individual [0 180] oscillation
    [~,b] = find(isnan(fp2ph{y,2})); % find first column with NaN's
    fp2ph{y,2}(:,[b(1):size(fp2ph{y,2},2)]) = []; % remove end of matrix with incomplete columns
    
    for z = 1:3
        fp2ph_state_norm{y,z} = normalize(fp2ph_state{y,z},1,'range');
        [~,b] = find(isnan(fp2ph_state{y,z})); % find first column with NaN's
        fp2ph_state_norm{y,z}(:,[b(1):size(fp2ph_state{y,z},2)]) = []; % remove end of matrix with incomplete columns
    end
    
end

%% DA fp norm to ACh phase
y = 1;
[~, edges] = histcounts(fp_deg(:,y), 'BinWidth', 10, 'BinLimits', [-180 180]); % genereate edges vector for phase, in degrees
[phIdx] = discretize(fp_deg(:,y), edges); % returns the indices of the bins that the elements of X fall into

y = 2;
da2achph_state = cell(3,1); for a = 1:3; da2achph_state{a} = nan(36, 3000); end % initialize matrix
da2achph_state_norm = cell(3,1); for a = 1:3; da2achph_state_norm{a} = nan(36, 3000); end % initialize matrix

for z = 1:length(unique(phIdx))
    ii = find(phIdx == z); % find index in signal that are this phase
    for c = 1:3
        ii_c = ii(ismember(ii, idx_cell{c})); % include index for this behavioral state
        da2achph_state{c}(z,[1:length(ii_c)]) = fp_mat(ii_c,y); % extract photometry values for this behavioral state
    end
end
for z = 1:3
    da2achph_state_norm{z} = normalize(da2achph_state{z},1,'range');
    [~,b] = find(isnan(da2achph_state_norm{z})); % find first column with NaN's
    da2achph_state_norm{z}(:,[b(1):size(da2achph_state{z},2)]) = []; % remove end of matrix with incomplete columns
end

%% PLOT EXAMPLE

fig = figure; fig.Position([3 4]) = [1375 800];
clr = {'r','g','b'};

for y = 1:2
    switch y; case 1; lbl = 'ACh'; case 2; lbl = 'DA'; end
    
    a = y; subplot(2,3,a); hold on
    for z = 1:3
        shadederrbar( edges(2:end)-5, nanmean(fp2ph_state{y,z},2), SEM(fp2ph_state{y,z},2), clr{z});
    end
    xlabel(sprintf('%s Phase',lbl)); xlim([-180 180]); xticks([-180:90:180]);
    ylabel(sprintf('%s Fluorescence (%dF/F)',lbl));
    legend({'imm','','mov','','rew'});
    title(sprintf('%s fp to %s phase',lbl,lbl)); axis('square');
    
    a = y+3; subplot(2,3,a); hold on
    for z = 1:3
        shadederrbar( edges(2:end)-5, nanmean(fp2ph_state_norm{y,z},2), SEM(fp2ph_state_norm{y,z},2), clr{z});
    end
    xlabel(sprintf('%s Phase',lbl)); xlim([-180 180]); xticks([-180:90:180]);
    ylabel(sprintf('%s Fluorescence (norm)',lbl)); ylim([0 1]);
    legend({'imm','','mov','','rew'});
    title(sprintf('%s fp norm to %s phase',lbl,lbl));axis('square');
end
subplot(2,3,6); hold on
for z = 1:3
    shadederrbar( edges(2:end)-5, nanmean(da2achph_state_norm{z},2), SEM(da2achph_state_norm{z},2), clr{z});
end
xlabel('ACh Phase'); xlim([-180 180]); xticks([-180:90:180]);
ylabel('DA Fluorescence (norm)'); ylim([0 1]);
legend({'imm','','mov','','rew'});
title('DA fp norm to ACh phase');axis('square');
movegui(gcf, 'center')
