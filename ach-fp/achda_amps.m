good_rew = [10:16,20,22:26,28:29,32:35,40:42,44:46]; % good rew
beh = modAChDA(good_rew);

x = 13;
fp_mat = [beh(x).FP{1}, beh(x).FP{2}]; 
fp_mat = fp_mat - nanmean(fp_mat);
Fs = 50;

if isfield(beh,'reward')
    rew = beh(x).reward(:)/Fs;
    lick = beh(x).lick(:)/Fs;
    lick_repeat = [diff(lick.*1000) > 50]; % Identify licks that are <50ms after previous lick
    lick_sub = lick; lick_sub(1) = [];
    lick = [lick(1); lick_sub(lick_repeat)];
    bin = 1/1000; window = [0 1];
    peth = getClusterPETH(lick, rew, bin, window); % PETH: lick aligned to reward in 1 ms bins
    cts = peth.cts{1}; % Lick counts in 1ms bins for each reward trial
    rew_lick0 = find(sum(cts,1) == 0); % Find reward index where total licks within window is 0
    rew_prewlick = [];
    rew([rew_lick0, rew_prewlick]) = nan; 
    rew = rew(~isnan(rew)).*Fs; % return to samples
    idx_rew = extractEventST([1:length(fp_mat)]', floor(rew), floor(rew)+49, 1); % identify sample during reward
else; idx_rew = []; end
idx_mov = extractEventST([1:length(fp_mat)]', beh(x).on, beh(x).off, 1); % identify sample during locomotion
idx_mov_nonRew = idx_mov(~ismember(idx_mov, idx_rew)); % exclude reward, include locomotion
idx_imm = extractEventST([1:length(fp_mat)]', beh(x).onRest, beh(x).offRest, 1); % identify sample during locomotion
idx_imm_nonRew = idx_imm(~ismember(idx_imm, idx_rew)); % exclude reward, include rest
idx_c = cell(3,1); idx_c{1} = idx_imm_nonRew; idx_c{2} = idx_mov_nonRew; idx_c{3} = idx_rew; % index into cell array for ease of iteration

%
binMat = cell(3,1); % chop data into 1s bins
for a = 1:3
    binNum = floor(length(idx_c{a})/Fs);
    sig = fp_mat(idx_c{a},1); % extract ACh signal
    sigBin = reshape(sig(1:Fs*binNum), [Fs, binNum]); % reshape data into 1s bins
    sigBinMin = min(sigBin); % minimum ACh in each 1s bin
    
    sig = fp_mat(idx_c{a},2); % extract DA signal
    sigBin = reshape(sig(1:Fs*binNum), [Fs, binNum]); % reshape data into 1s bins
    sigBinMax = max(sigBin); % maximum DA in each 1s bin
    
    binMat{a} = [sigBinMin(:), sigBinMax(:)];
end

%%
fig = figure; fig.Position(3) = 1375;
colors = [1 0 0 0.1; 0 1 0 0.1; 0 0 1 0.1]; 
edges_da = [floor(min(fp_mat(:,1))) : 0.5 : ceil(max(fp_mat(:,1)))];
edges_ach = [floor(min(fp_mat(:,2))) : 0.5 : ceil(max(fp_mat(:,2)))];

subplot(1,3,1); hold on
for a = [1 2 3]
    histogram(binMat{a}(:,1), 'BinEdges', edges_ach, 'Normalization', 'probability', 'FaceColor', colors(a,[1:3]), 'FaceAlpha', 0.2);
end
xlabel('ACh minima'); ylabel('probability'); title('ACh minimum');

subplot(1,3,2); hold on
for a = [1 2 3]
    histogram(binMat{a}(:,2), 'BinEdges', edges_da, 'Normalization', 'probability', 'FaceColor', colors(a,[1:3]), 'FaceAlpha', 0.2);
end
legend({'immobility','locomotion','reward'});
xlabel('DA maxima'); ylabel('probability'); title('DA maximum');

subplot(1,3,3); hold on
for a = 1:3
    plot(binMat{a}(:,1), binMat{a}(:,2), '.', 'color', colors(a,:), 'MarkerSize', 10);
end
title(sprintf('%s',beh(x).rec)); axis equal
xlabel('ACh minimum'); ylabel('DA maximum');

movegui(gcf,'center');