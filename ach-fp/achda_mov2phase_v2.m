%% summary
% Find instantaneous phase of each signal
% Plot distribution of %dF/F photometry signal w.r.t. instantaneous phase:
%   (A) for each signal to its own phase, comparing behavioral states
%   (B) for each signal, normalized, to its own phase to demonstrate
%   relationship with focus on temporal aspect
%   (C) for DA photometry signal, normalized, to the phase of ACh

%%
% beh = modAChDA; beh(42) = [];
win = [-2 2]; % window for alignment of signal to peaks

h = waitbar(0, 'movement to photometry phase peak');
mat = struct;
for x = 1:length(beh)
    %% 
    mov_mat = [beh(x).vel, getAcc(beh(x).vel)]; % extract velocity, acceleration
    nSamp = length(beh(x).time);
    fp_mat = []; fp_phase = []; fp_deg = []; % clear photometry 
    peakIdxOnset = []; peakIdxLag = []; % clear peaks
    mov2fpk_onset = cell(2,2); mov2fpk_mov = cell(2,2); mov2fpk_shuff = cell(2,2); % initialie output 
    
    %% Index of behavioral states
    idx_cell = cell(3,1); 
    if isfield(beh,'reward')
        idx_rew = extractEventST([1:nSamp]', floor(beh(x).reward), floor(beh(x).reward)+100, 1); % identify sample during reward
        idx_cell{3} = idx_rew; 
    else; idx_rew = []; idx_cell = cell(2,1); end
    idx_mov = extractEventST([1:nSamp]', beh(x).on, beh(x).off, 1); % identify sample during locomotion
    idx_mov_nonRew = idx_mov(~ismember(idx_mov, idx_rew)); % exclude reward, include locomotion
    idx_imm = extractEventST([1:nSamp]', beh(x).onRest, beh(x).offRest, 1); % identify sample during locomotion
    idx_imm_nonRew = idx_imm(~ismember(idx_imm, idx_rew)); % exclude reward, include rest
    idx_cell{1} = idx_imm_nonRew; idx_cell{2} = idx_mov_nonRew; % index into cell array for ease of iteration

    for y = 1:2
        fp_mat(:,y) = beh(x).FP{y} - nanmean(beh(x).FP{y}); 

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
        fp_phase = data_phase;
        fp_deg = rad2deg(data_phase);

        %Take the absolute of the filtered signal and calculate the standard deviation
        rmssig  = abs(data_filt);
        stdsig = std(rmssig);

        %Find the index of peaks that are more than NumStd standard deviations (NumStd could be 1.5 standard deviations)
        % finds the 0 degree phase indices that cross 1.5 standard deviations. Might have to remove the ; depending on whether your data is a row or column array
        NumStd = 1.5;
        peakIdx = find(data_filt>(NumStd*stdsig) & [0; diff(data_phase>0)]);      
        peakIdxMov = peakIdx(ismember(peakIdx, idx_mov_nonRew));
        
        %% First peak within a movement
        for z = 1:length(beh(x).on)
            [~,ii] = min(abs(peakIdx - beh(x).on(z))); % find index of closest peakIdx to movement onset
            if peakIdx(ii) - beh(x).on(z) < 0
                ii = ii + 1; % adjust if index is for peakIdx that is before onset
            end
            peakIdxOnset(z,y) = peakIdx(ii); % peakIdx of 1st peak within a movement bout
            peakIdxLag(z,y) = peakIdx(ii) - beh(x).on(z); % latency of peak, in samples, from onset
            if peakIdxLag(z,y) > 50
                peakIdxOnset(z,y) = nan; peakIdxLag(z,y) = nan;
            end
        end
        
        %% Shuffle peaks
        peakIdxShuff = shuffleST(peakIdx, 1); peakIdxShuff = peakIdxShuff{1};
        peakIdxShuffMov = peakIdxShuff(ismember(peakIdxShuff, idx_mov_nonRew));

        %% Align movmement to peak
        for z = 1:2
           [mov2fpk_onset{z,y}, sta_time] = getSTA(mov_mat(:,z), peakIdxOnset(:,y)/Fs, Fs, win); % align velocity or acceleration to 1st peak within movement bout
           [mov2fpk_mov{z,y}, sta_time] = getSTA(mov_mat(:,z), peakIdxMov/Fs, Fs, win); % align velocity or acceleration to 1st peak within movement bout
           mov2fpk_shuff{z,y} = getSTA(mov_mat(:,z), peakIdxShuffMov/Fs, Fs, win);
        end
        
    end
    %% Save into output structure
    mat(x).rec = beh(x).rec;
    mat(x).mov2fpk_mov = mov2fpk_mov; % {1,1} vel to ACh | {2,1} acc to ACh | {1,2} vel to DA | {2,2} acc to DA
    mat(x).mov2fpk_onset = mov2fpk_onset; % {1,1} vel to ACh | {2,1} acc to ACh | {1,2} vel to DA | {2,2} acc to DA
    mat(x).mov2fpk_shuff = mov2fpk_shuff;
    mat(x).peakIdxOnset = peakIdxOnset;
    mat(x).peakIdxLag = peakIdxLag;
    %%
    waitbar(x/length(beh),h);
end
fprintf('Done! \n'); close(h);

%% testing
x = 4; z = 1;
r = [beh(x).on(z):beh(x).off(z)]; 
r2 = mat(x).peakIdxOnset(z,:);
figure;
sp(1) = subplot(2,1,1); hold on
plot(beh(x).time(r),beh(x).FP{1}(r),'g');
plot(beh(x).time(r),beh(x).FP{2}(r),'m');
plot(beh(x).time(r2(1)),beh(x).FP{1}(r2(1)),'.g','MarkerSize',20);
plot(beh(x).time(r2(2)),beh(x).FP{2}(r2(2)),'.m','MarkerSize',20);
ylabel('Fluorescence');
title(sprintf('%s - bout#%d',beh(x).rec,z));

sp(2) = subplot(2,1,2); hold on
plot(beh(x).time(r),beh(x).vel(r),'k');
plot([beh(x).time(r2(1)) beh(x).time(r2(1))],[0 max(beh(x).vel(r))],'-g');
plot([beh(x).time(r2(2)) beh(x).time(r2(2))],[0 max(beh(x).vel(r))],'-m');
ylabel('Velocity (cm/s)');
yyaxis right; acc = getAcc(beh(x).vel); plot(beh(x).time(r),acc(r),'b'); ylabel('Acc (cm/s^2');

%% average across movements for each recording
% mov2fpk: {1,1} vel to ACh | {2,1} acc to ACh | {1,2} vel to DA | {2,2} acc to DA
mov2fpk_rec = cell(2,2); mov2fpk_shuff = cell(2,2);

for x = 1:length(mat)
    for a = 1:2; for b = 1:2
            mov2fpk_rec{a,b}(:,x) = nanmean(mat(x).mov2fpk_mov{a,b},2); 
            % base = nanmean(mov2fpk_rec{a,b}([find(sta_time == -2):find(sta_time == -1)],x));
            base = mov2fpk_rec{a,b}([find(sta_time == -1)],x);
            mov2fpk_rec{a,b}(:,x) = mov2fpk_rec{a,b}(:,x) - base;
            
            mov2fpk_shuff{a,b}(:,x) = nanmean(mat(x).mov2fpk_shuff{a,b},2);
            mov2fpk_shuff{a,b}(:,x) = mov2fpk_shuff{a,b}(:,x) - mov2fpk_shuff{a,b}([find(sta_time == -1)],x);
        end; end
end

% plot average across all recordings
figure;
sp_v(1) = subplot(2,2,1); hold on
a = 1; b = 1;
shadederrbar(sta_time, nanmean(mov2fpk_rec{a,b},2), SEM(mov2fpk_rec{a,b},2), 'g');
shadederrbar(sta_time, nanmean(mov2fpk_shuff{a,b},2), SEM(mov2fpk_shuff{a,b},2), 'k');
plot([0 0],[-1 1],'--k');
xlabel('latency from ACh peak (s)'); xlim([-1 1]); 
ylabel('change in velocity (cm/s)');
title('Velocity to ACh (n = 46 rec)'); axis('square');

sp_v(2) = subplot(2,2,2); hold on 
a = 1; b = 2;
shadederrbar(sta_time, nanmean(mov2fpk_rec{a,b},2), SEM(mov2fpk_rec{a,b},2), 'm');
shadederrbar(sta_time, nanmean(mov2fpk_shuff{a,b},2), SEM(mov2fpk_shuff{a,b},2), 'k');
plot([0 0],[-1 1],'--k');
xlabel('latency from DA peak (s)'); xlim([-1 1]); 
ylabel('change in velocity (cm/s)');
title('Velocity to DA '); axis('square');

sp_a(1) = subplot(2,2,3); hold on 
a = 2; b = 1;
shadederrbar(sta_time, nanmean(mov2fpk_rec{a,b},2), SEM(mov2fpk_rec{a,b},2), 'g');
shadederrbar(sta_time, nanmean(mov2fpk_shuff{a,b},2), SEM(mov2fpk_shuff{a,b},2), 'k');
plot([0 0],[-0.1 0.1],'--k');
xlabel('latency from ACh peak (s)'); xlim([-1 1]); 
ylabel('change in Acc (cm/s^2)');
title('Acceleration to ACh'); axis('square');

sp_a(2) = subplot(2,2,4); hold on 
a = 2; b = 2;
shadederrbar(sta_time, nanmean(mov2fpk_rec{a,b},2), SEM(mov2fpk_rec{a,b},2), 'm');
shadederrbar(sta_time, nanmean(mov2fpk_shuff{a,b},2), SEM(mov2fpk_shuff{a,b},2), 'k');
plot([0 0],[-0.1 0.1],'--k');
xlabel('latency from DA peak (s)'); xlim([-1 1]); 
ylabel('change in Acc (cm/s^2)');
title('Acceleration to DA'); axis('square');

linkaxes(sp_a,'y'); linkaxes(sp_v,'y');
%% average across movements for each recording
% mov2fpk: {1,1} vel to ACh | {2,1} acc to ACh | {1,2} vel to DA | {2,2} acc to DA
mov2fpk_rec = cell(2,2);

for x = 1:length(mat)
    for a = 1:2; for b = 1:2
            mov2fpk_rec{a,b}(:,x) = nanmean(mat(x).mov2fpk_onset{a,b},2); 
            base = nanmean(mov2fpk_rec{a,b}([find(sta_time == -4):find(sta_time == -1)],x));
            mov2fpk_rec{a,b}(:,x) = (mov2fpk_rec{a,b}(:,x) - base)./base;
            mov2fpk_rec{a,b}(:,x) = mov2fpk_rec{a,b}(:,x).*100;
        end; end
end

% plot average across all recordings
fig = figure; fig.Position(3) = 1000;
subplot(1,2,1); hold on
shadederrbar(sta_time, nanmean(mov2fpk_rec{1,1},2), SEM(mov2fpk_rec{1,1},2), 'g');
plot([0 0],[-10 30],'--k');
xlabel('latency from ACh peak (s)'); xlim([-1 1]);
ylabel('delta Velocity (%)');
title('Velocity to ACh (n = 46 rec)'); axis('square');

subplot(1,2,2); hold on 
shadederrbar(sta_time, nanmean(mov2fpk_rec{1,2},2), SEM(mov2fpk_rec{1,2},2), 'm');
plot([0 0],[-5 15],'--k');
xlabel('latency from DA peak (s)'); xlim([-1 1]);
ylabel('delta Velocity (%)');
title('Velocity to DA'); axis('square');