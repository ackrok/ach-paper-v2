y = 1; beh = s(y).s;

align_rew = plot_fp2event(beh,[-6 2],0); % reward alignment
% [align_on,t_on] = plot_fp2event(beh,[-6 2],0); % onset alignment

%%
prc_captured = [];
ach2rew = cell(length(beh),2);
ach2rew_acsf = cell(2,1);

for x = 1:length(beh)
 
    fp = beh(x).FP{1}; Fs = beh(x).Fs; % extract photometry signal from structure
    fp_mu = nanmean(fp); % mean of entire photometry signal
    fp = fp - fp_mu; % SUBTRACT baseline% from fp%, now centered on BASELINE
    
    %%
    rew_true = beh(x).reward(~isnan(align_rew{x,y}(1,:))); % excluded non-rewarded trials
    idx_rew = extractEventST([1:length(fp)]', floor(rew_true), floor(rew_true)+100, 1); % identify sample during reward
    idx_inf = [s(y).win(x,1).*(50*60), s(y).win(x,2).*(50*60)]';
    idx_rew(idx_rew < idx_inf(1)) = []; idx_rew(idx_rew > idx_inf(2)) = []; % exclude rewards outside of infusion window
    rew_true_inf = rew_true; rew_true_inf(rew_true_inf < idx_inf(1)) = []; rew_true_inf(rew_true_inf > idx_inf(2)) = [];
    
    thres = s(y).thres(x,1); width = s(y).thres(x,1); % THRES/WIDTH for pause
    fp_rew = fp(idx_rew); % photometry during reward window
    [rew_idxMax, rew_idxGood, rew_valMaxAll, rew_crossGroups, rew_crossStart] = findIdxPause(fp_rew, 'pause', thres, width); 
    rew_startGood = rew_crossStart(rew_idxGood,:); % Extract start/end times for good pauses
    rew_valGood = rew_valMaxAll(rew_idxGood);
    tmp = [1:101:length(fp_rew)]'; % Reward start index, in index values of fp_rew
    tmp = ceil(rew_startGood./101); % Divide by length of window to find which reward this pause is in
    tmp((diff(tmp,1,2) ~= 0),:) = nan; % NaN if pause start/end not within same reward window
    uni = unique(tmp); uni = uni(~isnan(uni)); % Unique reward windows that have crossings
    rew_valGoodUnique = nan(length(uni),1); % Initialize
    rew_startGoodUnique = nan(length(uni),2); % Initialize
    rew_idxMaxUnique = nan(length(uni),1); % Initialize
    for c = 1:length(uni)
        cc = find(tmp(:,1) == uni(c)); % idxGood that matches unique reward windows
        [rew_valGoodUnique(c),ii] = min(rew_valGood(cc)); % Keep minimum for multiple crosses in single reward period
        rew_startGoodUnique(c,:) = rew_startGood(cc(ii),:); 
        rew_idxMaxUnique(c) = rew_idxMax(cc(ii));
    end
    rew_valMax{x} = rew_valGoodUnique;
    prc_captured(x,1) = length(rew_idxMaxUnique)/length(rew_true_inf);
    
    [a,sta_time] = getSTA(fp(idx_rew), rew_idxMaxUnique/Fs, Fs, [-4 2]);
    s(y).ach2rew{x,1} = a;
    ach2rew{x,1} = a;
    a = a - nanmean(a(find(sta_time == -4):find(sta_time == -1),:)); % subtract baseline [-4 -1]
    a = a(find(sta_time == -1):find(sta_time == 1),:); % restrict to [-1 1] window
    a = nanmean(a,2);
    b = (a - min(a))/(a(1,:) - min(a)); % min max normalization
    ach2rew_acsf{1}(:,x) = b;
    
    %%
    idx_onset = extractEventST([1:length(fp)]', beh(x).on, beh(x).on+25, 1); % Determine indices for onset window: [0 +0.5]
    idx_onset(idx_onset < idx_inf(1)) = []; idx_onset(idx_onset > idx_inf(2)) = []; % exclude onsets outside of infusion window
    on_true = beh(x).on;
    on_true_inf = on_true; on_true_inf(on_true_inf < idx_inf(1)) = []; on_true_inf(on_true_inf > idx_inf(2)) = [];
    idx_true_inf = find(ismember(beh(x).on, on_true_inf));
    
    thres = s(y).thres(x,2); 
    width = s(y).thres(x,2); % THRES/WIDTH for peak
    ev = on_true_inf./Fs; % For only trials during infusion window
    a_on = getSTA(fp, ev, Fs, [0 0.5]);
    base = nanmean(getSTA(fp, ev, Fs, [-4 -2]),1); % Subtract baseline
    a_on = a_on - base;
    a_on = a_on - a_on(1,:);
    a_on = a_on(:,all(~isnan(a_on))); % remove NaN trials
    [m, i] = max(a_on); % i = i+find(t_on == 0)-1; % onset window: [0 1]s
    cross = find(m > thres); % Find onset trials where onset peak exceeds threshold
    peak_idx = on_true_inf(cross) + i(cross)' - 1;
    peak_val = m(cross)';
    
    a = getSTA(fp, peak_idx./Fs, Fs, [-4 2]);
    ach2rew{x,2} = a;
    s(y).ach2rew{x,2} = a;
    a = a - nanmean(a(find(sta_time == -4):find(sta_time == -1),:)); % subtract baseline [-4 -1]
    a = a(find(sta_time == -1):find(sta_time == 1),:); % restrict to [-1 1] window
    a = nanmean(a,2);
    b = (a - a(1,:))/(a(51,:) - a(1,:)); % min max normalization
    ach2rew_acsf{2}(:,x) = b;
end