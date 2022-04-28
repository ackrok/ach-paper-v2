%% summary
% Find instantaneous phase of each signal
% Plot distribution of %dF/F photometry signal w.r.t. instantaneous phase:
%   (A) for each signal to its own phase, comparing behavioral states
%   (B) for each signal, normalized, to its own phase to demonstrate
%   relationship with focus on temporal aspect
%   (C) for DA photometry signal, normalized, to the phase of ACh

%%
beh = modAChDA; beh([27 33 42]) = [];

vel2ph_beh = cell(2,3); 
acc2ph_beh = cell(2,3);

h = waitbar(0, 'photometry to instantaneous phase');
for x = 1:length(beh)

    fp_mat = []; fp_phase = []; fp_deg = [];
    vel_mat = [beh(x).vel]; % vel_mat = (vel_mat - mode(vel_mat))./(max(vel_mat) - mode(vel_mat));
    acc_mat = [getAcc(beh(x).vel)];
    
    mov2ph = cell(2,2); for a = 1:2; for b = 1:2; mov2ph{a,b} = nan(36, 3000); end; end % initialize matrix
    vel2ph_state = cell(2,3); for a = 1:2; for b = 1:3; vel2ph_state{a,b} = nan(36, 3000); end; end % initialize matrix
    acc2ph_state = vel2ph_state; % initialize matrix
    
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

        %% Find photometry "peaks", when inst. phase = 0
        %Bandpass filter
        signal = fp_mat(:,y);
        Fpass = [0.5 4];
        Fs = 50; %sampling rate, has to be at least double of your high pass frequency
        Wn = [Fpass(1)/(Fs/2) Fpass(2)/(Fs/2)]; 
        [b,a] = butter(3,Wn);
        data_filt= filtfilt(b,a,signal); % signal is your lfp data, output is the filtered data

        %Instantaneous phase
        H = hilbert(double(data_filt));
        data_phase = angle(H); % output is the instantaneous phase
        fp_phase(:,y) = data_phase;
        fp_deg(:,y) = rad2deg(data_phase);

        %Take the absolute of the filtered signal and calculate the standard deviation
        rmssig  = abs(data_filt);
        stdsig = std(rmssig);

        %Find the index of peaks that are more than NumStd standard deviations (NumStd could be 1.5 standard deviations)
        % finds the 0 degree phase indices that cross 1.5 standard deviations. Might have to remove the ; depending on whether your data is a row or column array
        NumStd = 2;
        peakIdx = find(data_filt>(NumStd*stdsig) & [0; diff(data_phase>0)]);         

        %% STA: signal to peaks
        win = [-1 1];
        [mov2ph{y,1},sta_time] = getSTA(vel_mat, peakIdx./Fs, Fs, win);
        mov2ph{y,2} = getSTA(acc_mat, peakIdx./Fs, Fs, win);
        for z = 1:3
            ii_c = peakIdx(ismember(peakIdx, idx_cell{z})); % include index for this behavioral state
            vel2ph_state{y,z} = getSTA(vel_mat, ii_c./Fs, Fs, win);
            acc2ph_state{y,z} = getSTA(acc_mat, ii_c./Fs, Fs, win);
            vel2ph_beh{y,z}(:,x) = nanmean(vel2ph_state{y,z},2);
            acc2ph_beh{y,z}(:,x) = nanmean(acc2ph_state{y,z},2);
        end
        
        %% Align signal to phase
%         [~, edges] = histcounts(fp_deg(:,y), 'BinWidth', 10, 'BinLimits', [-180 180]); % genereate edges vector for phase, in degrees
%         [phIdx] = discretize(fp_deg(:,y), edges); % returns the indices of the bins that the elements of X fall into
% 
%         for z = 1:3
%             for p = 1:length(unique(phIdx))
%                 ii = find(phIdx == p); % find index in signal that are this phase
%                 mov2ph{y,1}(p,[1:length(ii)]) = vel_mat(ii); % extract signal values matching this phase
%                 mov2ph{y,2}(p,[1:length(ii)]) = acc_mat(ii); % extract signal values matching this phase
%                 
%                 ii_c = ii(ismember(ii, idx_cell{z})); % include index for this behavioral state
%                 vel2ph_state{y,z}(p,[1:length(ii_c)]) = vel_mat(ii_c); % extract values for this behavioral state
%                 acc2ph_state{y,z}(p,[1:length(ii_c)]) = acc_mat(ii_c); % extract values for this behavioral state
%             end
%             % remove any columns with NaNs:
%             % [~,b] = find(isnan(mov2ph{y,z})); % find first column with NaN's
%             % mov2ph{y,z}(:,[b(1):size(mov2ph{y,z},2)]) = []; % remove end of matrix with incomplete columns
%             [~,b] = find(isnan(vel2ph_state{y,z})); % find first column with NaN's
%             vel2ph_state{y,z}(:,[b(1):size(vel2ph_state{y,z},2)]) = []; % remove end of matrix with incomplete columns
%             [~,b] = find(isnan(acc2ph_state{y,z})); % find first column with NaN's
%             acc2ph_state{y,z}(:,[b(1):size(acc2ph_state{y,z},2)]) = []; % remove end of matrix with incomplete columns
%             
%             % save into output cell array:
%             vel2ph_beh{y,z}(:,x) = nanmean(vel2ph_state{y,z},2);
%             acc2ph_beh{y,z}(:,x) = nanmean(acc2ph_state{y,z},2);
%         end
        
    end

    %%
    waitbar(x/length(beh),h);
end
fprintf('Done! \n'); close(h);
% clearvars -except modAChDA beh edges fp2ph_beh fp2ph_norm_beh da2achph_norm_beh

%% N = X mice: DA to ACh pause
tmp = {}; for x = 1:length(beh); tmp{x} = strtok(beh(x).rec,'-'); end
uni = unique(tmp); nAn = length(uni);
vel2ph_an = cell(2,3); % Initialize cell arrays
acc2ph_an = cell(2,3); 

for x = 1:nAn
    ii = find(strcmp(tmp,uni{x}));
    for z = 1:3
        for y = 1:2
            vel2ph_an{y,z}(:,x) = nanmean([vel2ph_beh{y,z}(:,ii)],2);
            % vel2ph_an{y,z}(:,x) = normalize(nanmean([vel2ph_beh{y,z}(:,ii)],2),'range');
            acc2ph_an{y,z}(:,x) = nanmean([acc2ph_beh{y,z}(:,ii)],2);
            % acc2ph_an{y,z}(:,x) = normalize(nanmean([acc2ph_beh{y,z}(:,ii)],2),'range');
        end
    end
end

%% PLOT by ANIMAL
fig = figure; fig.Position([3 4]) = [1000 800]; movegui(gcf, 'center')
clr = {'r','g','b'};
mid = edges(2:end)-5;

for y = 1:2
    switch y; case 1; lbl = 'ACh'; case 2; lbl = 'DA'; end
    
    a = y; sp(y) = subplot(2,2,a); hold on
    for z = 1:2
        shadederrbar( sta_time, nanmean(vel2ph_an{y,z},2), SEM(vel2ph_an{y,z},2), clr{z});
    end
    xlabel(sprintf('lag to %s peak (s)',lbl)); xlim([-1 1]);
    ylabel('Velocity');
    legend({'mov'});
    title(sprintf('Velocity to %s phase (n = %d)',lbl,size(vel2ph_an{1,1},2))); axis('square');
    
    a = y+2; sp2(y) = subplot(2,2,a); hold on
    for z = 1:2
        shadederrbar( sta_time, nanmean(acc2ph_an{y,z},2), SEM(acc2ph_an{y,z},2), clr{z});
    end
    plot([0 0],[-0.1 0.1],'--k');
    xlabel(sprintf('lag to %s peak (s)',lbl)); xlim([-1 1]);
    ylabel('Acceleration');
    legend({'imm','','mov'});
    title(sprintf('Acceleration to %s phase',lbl));axis('square');
end
linkaxes(sp,'y'); linkaxes(sp2,'y');
% subplot(2,3,6); hold on
% for z = 1:3
%     % shadederrbar( edges(2:end)-5, nanmean(da2achph_norm_an{1,z},2), SEM(da2achph_norm_an{1,z},2), clr{z});
%     plot( edges(2:end)-5, nanmean(da2achph_norm_an{1,z},2), clr{z});
% end
% plot( edges(2:end)-5, nanmean(fp2ph_norm_an{1,1},2), '--k');
% xlabel('ACh Phase'); xlim([-180 180]); xticks([-180:90:180]);
% ylabel('DA Fluorescence (norm)'); ylim([0 1]);
% % legend({'imm','','mov','','rew','','ACh'});
% title('DA fp norm to ACh phase');axis('square');
% movegui(gcf, 'center')
