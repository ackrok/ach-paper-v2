%% summary
% Find instantaneous phase of each signal
% Plot distribution of %dF/F photometry signal w.r.t. instantaneous phase:
%   for ACh photometry signal, normalized, to the phase of LFP

%% (1) Load data
load('C:\Users\Anya\Desktop\IV_LOCAL\data_comb\fullbeh_ACh_v2')
beh = behACh; % Extract ACh recordings
load('C:\Users\Anya\Desktop\IV_LOCAL\data_comb\fullcin_Feb21_ACh+DA+PF');
sub = cinwt;  % Extract CINs from recordings with ACh

%% (2) PETH
lbl = 'peak'; % CHANGE

Fs = 50; nShuff = 5;
bin = 0.02; window = [-1 1]; %CHANGE: window for PETH

mat = struct; % Initialize structure
h = waitbar(0, 'PETH: CIN spikes to FP peak times'); fprintf('Running ... ');
for x = 1:length(beh)
    fprintf('%s ... ',beh(x).rec);

    %% Extract spike times
    idx = find(strcmp({sub.rec},beh(x).rec));
    if isempty(idx); continue; end
    st = {sub(idx).st}; % Extract spike times of units from this recording
    fr = []; % Compute unit firing rate during movement and rest
    for y = 1:length(idx)
        % fr(y,1) = 1/mean(diff(extractEventST(st{y},beh(x).on./Fs,beh(x).off./Fs,0))); % Event times during movement
        fr(y,2) = 1/mean(diff(extractEventST(st{y},beh(x).onRest./Fs,beh(x).offRest./Fs,0))); % Event times during rest
        % fr(y,3) = 1/mean(diff(st{y}));
    end
    
    %% Compute fp peaks
    fp = beh(x).FP{1} - nanmean(beh(x).FP{1}); % PHOTOMETRY

    %% Identify fp peaks
    % Bandpass filter
    signal = fp;
    Fpass = [0.5 4];
    Fs = 50; %sampling rate, has to be at least double of your high pass frequency
    Wn = [Fpass(1)/(Fs/2) Fpass(2)/(Fs/2)]; 
    [b,a] = butter(3,Wn);
    data_filt = filtfilt(b,a,signal); % signal is your lfp data, output is the filtered data
    H = hilbert(double(data_filt));
    data_phase = angle(H); % output is the instantaneous phase
    data_deg = rad2deg(data_phase);
    % Take the absolute of the filtered signal and calculate the standard deviation
    rmssig  = abs(data_filt);
    stdsig = std(rmssig);
    % Find the index of peaks that are more than NumStd standard deviations (NumStd could be 1.5 standard deviations)
    % finds the 0 degree phase indices that cross 1.5 standard deviations. Might have to remove the ; depending on whether your data is a row or column array
    NumStd = 1.5;
    idxPeak = find(data_filt>(NumStd*stdsig) & [0; diff(data_phase>0)]);         
    idxPause = find(-data_filt>(NumStd*stdsig) & [0; diff(-data_phase>0)]);
    
    switch lbl; case 'pause'; locs = idxPause./Fs; % convert to seconds
        case 'peak'; locs = idxPeak./Fs; end
    
    %% Extract event times during rest & movement
    ev_sub = cell(1,3); %ev_sub_forShuff = cell(1,2);
    % ev_sub{1} = extractEventST(locs, beh(x).on./Fs, beh(x).off./Fs, 1); % Event times during movement
    ev_sub{2} = extractEventST(locs, beh(x).onRest./Fs, beh(x).offRest./Fs, 1); % Event times during rest
    % ev_sub{3} = locs;
%     ev_sub_forShuff{1} = extractEventST(locs, beh(x).on/Fs, beh(x).off/Fs, 1); % Event times during movement
%     ev_sub_forShuff{2} = extractEventST(locs, beh(x).onRest/Fs, beh(x).offRest/Fs, 1); % Event times during rest

        %% PETH
    align = cell(1,length(ev_sub)); alignDelta = cell(1,length(ev_sub)); 
    % alignZ = cell(1,length(ev_sub)); 
    shuff95 = align; shuff50 = align; shuff5 = align;
    for z = 1:length(ev_sub)
        if isempty(ev_sub{z}); continue; end
        peth = getClusterPETH(st, ev_sub{z}, bin, window); % PETH: spike times aligned to fp peaks
        align{z} = peth.fr;
        h2 = waitbar(0, 'PETH: st to fp locs');
        for y = 1:length(idx)
            alignDelta{z} = [alignDelta{z}, (peth.fr(:,y)-fr(y,z))./fr(y,z)]; % Delta firing rate change
            stShuff = shuffleST(st{y}, nShuff);
            peth_shuff = getClusterPETH(stShuff, ev_sub{z}, bin, window); %PETH: shuffled spike times aligned to fp peaks
            % mu = nanmean(nanmean(peth_shuff.fr,2)); sigma = nanmean(nanstd(peth_shuff.fr,[], 2)); 
            % tmp_z = (peth.fr(:,y) - mu)./sigma; % z-score signal using shuffled mu, sigma
            % alignZ{z} = [alignZ{z}, tmp_z]; % add to matrix
            prc = prctile(peth_shuff.fr,[5 50 95],2); %5th, 50th, 95th percentile of shuffled PETH
            shuff50{z} = [shuff50{z}, (prc(:,2)-fr(y,z))./fr(y,z)]; 
            shuff95{z} = [shuff95{z}, (prc(:,3)-fr(y,z))./fr(y,z)]; 
            shuff5{z} = [shuff5{z}, (prc(:,1)-fr(y,z))./fr(y,z)]; 
            waitbar(y/length(idx),h2);
        end
        close(h2);
    end
    
    %% STA
    sigal = fp;
    sta_fp = cell(1,length(ev_sub));
    for z = 1:length(ev_sub)
        if isempty(ev_sub{z}); continue; end
        [sta_fp{z},sta_time] = getSTA(sigal, ev_sub{z}, Fs, [window(1), window(2)]); % Align photometry to FP peaks
    end
    
    %% Load into output structure
    mat(x).rec = beh(x).rec; mat(x).FPnames{1} = beh(x).FPnames{1}; mat(x).n = [sub(idx).n]; mat(x).fr = fr;
    % mat(x).align_mvmt = align{1}; mat(x).alignDelta_mvmt = alignDelta{1}; 
    mat(x).align_rest = align{2}; mat(x).alignDelta_rest = alignDelta{2}; 
    % mat(x).align_full = align{3}; mat(x).alignDelta_full = alignDelta{3};
    % mat(x).alignZ_mvmt = alignZ{1}; mat(x).alignZ_rest = alignZ{2}; mat(x).alignZ_full = alignZ{3};
    % mat(x).shuff5_mvmt = shuff5{1}; mat(x).shuff50_mvmt = shuff50{1}; mat(x).shuff95_mvmt = shuff95{1};
    mat(x).shuff5_rest = shuff5{2}; mat(x).shuff50_rest = shuff50{2}; mat(x).shuff95_rest = shuff95{2}; 
    % mat(x).shuff5_full = shuff5{3}; mat(x).shuff50_full = shuff50{3}; mat(x).shuff95_full = shuff95{3};
    % mat(x).sta_mvmt = sta_fp{1}; 
    mat(x).sta_rest = sta_fp{2}; 
    % mat(x).sta_full = sta_fp{3};
   %%
    waitbar(x/length(beh),h); 
    fprintf('done. \n');
end
close(h); fprintf('\n Done: aligning spikes to ACh fp events \n');
time = peth.time; 

%% (3) Extract from output structure
% align_mvmt = []; alignDelta_mvmt = []; alignDelta_full = [];
align_rest = []; alignDelta_rest = []; 
%alignZ_mvmt = []; alignZ_rest = []; 
shuff50_rest = []; shuff95_rest = []; shuff5_rest = [];
for x = 1:length(mat)
    if isempty(mat(x).align_rest); continue; end
    % align_mvmt = [align_mvmt, mat(x).align_mvmt]; % Concatenate all units
    align_rest = [align_rest, mat(x).align_rest];
    % alignDelta_mvmt = [alignDelta_mvmt, mat(x).alignDelta_mvmt]; 
    alignDelta_rest = [alignDelta_rest, mat(x).alignDelta_rest];
    % alignDelta_full = [alignDelta_full, mat(x).alignDelta_full];
%     alignZ_mvmt = [alignZ_mvmt, mat(x).alignZ_mvmt];
%     alignZ_rest = [alignZ_rest, mat(x).alignZ_rest];
    shuff95_rest = [shuff95_rest, mat(x).shuff95_rest];
    shuff50_rest = [shuff50_rest, mat(x).shuff50_rest];
    shuff5_rest = [shuff5_rest, mat(x).shuff5_rest];
end

fprintf('DONE \n');