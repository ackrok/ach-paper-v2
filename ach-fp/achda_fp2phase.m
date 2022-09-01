%% summary
% Find instantaneous phase of each signal
% Plot distribution of %dF/F photometry signal w.r.t. instantaneous phase:
%   (A) for each signal to its own phase, comparing behavioral states
%   (B) for each signal, normalized, to its own phase to demonstrate
%   relationship with focus on temporal aspect
%   (C) for DA photometry signal, normalized, to the phase of ACh

%%
% beh = modAChDA; beh(42) = [];
NumStd = 1.5;

fp2ph_beh = cell(2,3);
fp2ph_norm_beh = fp2ph_beh;
da2achph_beh = fp2ph_beh;
da2achph_norm_beh = fp2ph_beh;

h = waitbar(0, 'photometry to instantaneous phase');
for x = 1:length(beh)

    fp_mat = []; fp_phase = []; fp_deg = [];
    fp2ph = cell(2,2); for a = 1:2; for b = 1:2; fp2ph{a,b} = nan(36, 10000); end; end % initialize matrix
    fp2ph_state = cell(2,3); for a = 1:2; for b = 1:3; fp2ph_state{a,b} = nan(36, 10000); end; end % initialize matrix
    fp2ph_state_norm = fp2ph_state; % initialize matrix

    for y = 1:2
        fp_mat(:,y) = beh(x).FP{y} - nanmean(beh(x).FP{y}); 

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
        % NumStd = 2; % MOVED ABOVE OUT OF FOR LOOP
        peakIdx = find(data_filt>(NumStd*stdsig) & [0; diff(data_phase>0)]);         

        %% Align photometry to phase
        [~, edges] = histcounts(fp_deg(:,y), 'BinWidth', 10, 'BinLimits', [-180 180]); % genereate edges vector for phase, in degrees
        [phIdx] = discretize(fp_deg(:,y), edges); % returns the indices of the bins that the elements of X fall into

        for z = 1:length(unique(phIdx))
            ii = find(phIdx == z); % find index in signal that are this phase
            fp2ph{y,1}(z,[1:length(ii)]) = fp_mat(ii,y); % extract photometry %dF/F values matching this phase
            for c = 1:length(idx_cell)
                ii_c = ii(ismember(ii, idx_cell{c})); % include index for this behavioral state
                fp2ph_state{y,c}(z,[1:length(ii_c)]) = fp_mat(ii_c,y); % extract photometry values for this behavioral state
            end
        end
        [~,b] = find(isnan(fp2ph{y,1})); % find first column with NaN's
        fp2ph{y,1}(:,[b(1):size(fp2ph{y,1},2)]) = []; % remove end of matrix with incomplete columns
        fp2ph{y,2} = normalize(fp2ph{y,1},1,'range'); % normalize over each individual [0 180] oscillation
        % [~,b] = find(isnan(fp2ph{y,2})); % find first column with NaN's
        % fp2ph{y,2}(:,[b(1):size(fp2ph{y,2},2)]) = []; % remove end of matrix with incomplete columns

        for z = 1:length(idx_cell)
            fp2ph_state_norm{y,z} = normalize(fp2ph_state{y,z},1,'range');
            [~,b] = find(isnan(fp2ph_state{y,z})); % find first column with NaN's
            fp2ph_state_norm{y,z}(:,[b(1):size(fp2ph_state{y,z},2)]) = []; % remove end of matrix with incomplete columns
        end

        %% save into output cell array
        for z = 1:length(idx_cell)
            fp2ph_beh{y,z}(:,x) = nanmean(fp2ph_state{y,z},2);
            fp2ph_norm_beh{y,z}(:,x) = nanmean(fp2ph_state_norm{y,z},2);
        end
    end

    %% DA fp norm to ACh phase
    for xx = 1:2
        switch xx 
            case 1
                y = 1;
                [~, edges] = histcounts(fp_deg(:,y), 'BinWidth', 10, 'BinLimits', [-180 180]); % genereate edges vector for phase, in degrees
                [phIdx] = discretize(fp_deg(:,y), edges); % returns the indices of the bins that the elements of X fall into
                y = 2;  
            case 2
                y = 2;
                [~, edges] = histcounts(fp_deg(:,y), 'BinWidth', 10, 'BinLimits', [-180 180]); % genereate edges vector for phase, in degrees
                [phIdx] = discretize(fp_deg(:,y), edges); % returns the indices of the bins that the elements of X fall into
                y = 1;      
        end
        da2achph_state = cell(3,1); for a = 1:3; da2achph_state{a} = nan(36, 10000); end % initialize matrix
        da2achph_state_norm = cell(3,1); for a = 1:3; da2achph_state_norm{a} = nan(36, 10000); end % initialize matrix
        for z = 1:length(unique(phIdx))
            ii = find(phIdx == z); % find index in signal that are this phase
            for c = 1:length(idx_cell)
                ii_c = ii(ismember(ii, idx_cell{c})); % include index for this behavioral state
                da2achph_state{c}(z,[1:length(ii_c)]) = fp_mat(ii_c,y); % extract photometry values for this behavioral state
            end
        end
        for z = 1:length(idx_cell)
            da2achph_state_norm{z} = normalize(da2achph_state{z},1,'range');
            [~,b] = find(isnan(da2achph_state_norm{z})); % find first column with NaN's
            da2achph_state_norm{z}(:,[b(1):size(da2achph_state{z},2)]) = []; % remove end of matrix with incomplete columns
            da2achph_norm_beh{xx,z}(:,x) = nanmean(da2achph_state_norm{z},2); % save into output cell array
            da2achph_beh{xx,z}(:,x) = nanmean(da2achph_state{z},2);
        end
    end

    %%
    waitbar(x/length(beh),h);
end
fprintf('Done! \n'); close(h);
% clearvars -except modAChDA beh edges fp2ph_beh fp2ph_norm_beh da2achph_norm_beh

% N = X mice: DA to ACh pause
tmp = {}; for x = 1:length(beh); tmp{x} = strtok(beh(x).rec,'-'); end
uni = unique(tmp); nAn = length(uni);
fp2ph_an = cell(2,3); fp2ph_norm_an = fp2ph_an; % Initialize cell arrays
da2achph_norm_an = fp2ph_an; da2achph_an = fp2ph_an; 

for x = 1:nAn
    ii = find(strcmp(tmp,uni{x}));
    for z = 1:length(idx_cell)
        for y = 1:2
            fp2ph_an{y,z}(:,x) = nanmean([fp2ph_beh{y,z}(:,ii)],2);
            fp2ph_norm_an{y,z}(:,x) = nanmean([fp2ph_norm_beh{y,z}(:,ii)],2);
            fp2ph_norm_an{y,z}(:,x) = normalize(fp2ph_norm_an{y,z}(:,x), 'range');
            da2achph_norm_an{y,z}(:,x) = nanmean([da2achph_norm_beh{y,z}(:,ii)],2);
            da2achph_norm_an{y,z}(:,x) = normalize(da2achph_norm_an{y,z}(:,x), 'range');
            da2achph_an{y,z}(:,x) = nanmean([da2achph_beh{y,z}(:,ii)],2);
        end
    end
end

%% PLOT by ANIMAL
% fig = figure; fig.Position([3 4]) = [1375 800];
% clr = {'r','g','b'};
% mid = edges(2:end)-5;
% for y = 1:2
%     switch y; case 1; lbl = 'ACh'; case 2; lbl = 'DA'; end
%     
%     a = y; sp(y) = subplot(2,3,a); hold on
%     for z = 1:length(idx_cell)
%         shadederrbar( mid, nanmean(fp2ph_an{y,z},2), SEM(fp2ph_an{y,z},2), clr{z});
%     end
%     xlabel(sprintf('%s Phase',lbl)); xlim([-180 180]); xticks([-180:90:180]);
%     ylabel(sprintf('%s Fluorescence (dF/F)',lbl));
%     legend({'imm','','mov','','rew'});
%     title(sprintf('%s fp to %s phase (n = %d)',lbl,lbl,size(fp2ph_an{1,1},2))); axis('square');
%     
%     a = y+2; subplot(2,3,a); hold on
%     for z = 1:length(idx_cell)
%         % shadederrbar( mid, nanmean(fp2ph_norm_an{y,z},2), SEM(fp2ph_norm_an{y,z},2), clr{z});
%         plot( mid, nanmean(fp2ph_norm_an{y,z},2), clr{z});
%     end
%     xlabel(sprintf('%s Phase',lbl)); xlim([-180 180]); xticks([-180:90:180]);
%     ylabel(sprintf('%s Fluorescence (norm)',lbl)); ylim([0 1]);
%     % legend({'imm','','mov','','rew'});
%     title(sprintf('%s fp norm to %s phase',lbl,lbl)); axis('square');
%     
%     a = y+4; subplot(2,3,a); hold on
%     sm = 5;
%     for z = 1:length(idx_cell)
%         shadederrbar( mid, movmean(nanmean(da2achph_norm_an{y,z},2),sm), movmean(SEM(da2achph_norm_an{y,z},2),sm), clr{z});
%         % plot( mid, movmean(nanmean(da2achph_norm_an{y,z},2),sm), clr{z});
%     end
%     switch y; case 1; xx = 2; case 2; xx = 1; end
%     plot( edges(2:end)-5, nanmean(fp2ph_norm_an{xx,1},2), '--k');
%     xlabel(sprintf('%s Phase',lbl)); xlim([-180 180]); xticks([-180:90:180]);
%     ylabel('Fluorescence (norm)'); ylim([0 1]);
%     % legend({'imm','','mov','','rew','','ACh'});
%     title(sprintf('fp norm to %s phase',lbl));axis('square');
% end
% 
% movegui(gcf, 'center')

%% QUANTIFICATION
% pref = cell(2,1);
% for a = 1:size(da2achph_norm_an,1); for b = 1:size(da2achph_norm_an,2)
%     [~,c] = max(da2achph_norm_an{a,b}); % maximum, normalized to 1 so this is the preferred phase
%     pref{a}(:,b) = mid(c)+180; % preferred phase at maximum
%     end; end
% 
% fig = figure; fig.Position(3) = 1000; clr = {'r','g','b'};
% for a = 1:2
%     subplot(1,2,a,polaraxes); hold on
%     for b = 1:3
%     polarhistogram(deg2rad(pref{a}(:,b)),6,'FaceColor',clr{b},'FaceAlpha',0.3,'Normalization','count');
%     end
% end