% good_rew = [10:16,20,22:26,28:29,32:35,40:42,44:46]; % good rew
good_rew = [10:14,20,22:26,28:29,32:35,40:42,44:46]; % good rew -- exclude later GZ03 trials
beh = modAChDA(good_rew);
%
bin = 0.5;
edges = [-15 : bin : 45]; % edges of full spectrum of photometry values
mat = struct;
for x = 1:length(beh)
    %%
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
        % sigBinMin = min(sigBin); % minimum ACh in each 1s bin
        sigBinMin = max(sigBin); % minimum ACh in each 1s bin
        
        sig = fp_mat(idx_c{a},2); % extract DA signal
        sigBin = reshape(sig(1:Fs*binNum), [Fs, binNum]); % reshape data into 1s bins
        sigBinMax = max(sigBin); % maximum DA in each 1s bin

        binMat{a} = [sigBinMin(:), sigBinMax(:)]; % save ACh minima, DA maxima
    end

    mat(x).rec = beh(x).rec; % add to output structure
    mat(x).binMat = binMat; % add to output structure

    %% PLOT EXAMPLE
    %{
    fig = figure; fig.Position(3) = 1375;
    colors = [1 0 0 0.1; 0 1 0 0.1; 0 0 1 0.1];
    y = 1; 
    subplot(1,3,y); hold on
    vec = [];
    for a = [1 2 3]
        histogram(binMat{a}(:,y), 'BinEdges', edges, 'Normalization', 'probability', 'FaceColor', colors(a,[1:3]), 'FaceAlpha', 0.2, 'EdgeAlpha', 0);
        n = histcounts(binMat{a}(:,y), edges, 'Normalization', 'probability');
        vec = [vec, n(:)];
    end
    [~,p] = kstest2(vec(:,1),vec(:,3)); [~,p(2)] = kstest2(vec(:,2),vec(:,3));
    xlabel('ACh minima'); ylabel('probability'); 
    title(sprintf('ACh minimum KS p = %1.4f, p = %1.4f',p(1),p(2))); axis square;

    y = 2; 
    subplot(1,3,y); hold on
    vec = [];
    for a = [1 2 3]
        histogram(binMat{a}(:,y), 'BinEdges', edges, 'Normalization', 'probability', 'FaceColor', colors(a,[1:3]), 'FaceAlpha', 0.2, 'EdgeAlpha', 0);
        n = histcounts(binMat{a}(:,y), edges, 'Normalization', 'probability');
        vec = [vec, n(:)];
    end
    [~,p] = kstest2(vec(:,1),vec(:,3)); [~,p(2)] = kstest2(vec(:,2),vec(:,3));
    legend({'immobility','locomotion','reward'});
    xlabel('DA maxima'); ylabel('probability'); 
    title(sprintf('DA maximum KS p = %1.4f, p = %1.4f',p(1),p(2))); axis square;

    subplot(1,3,3); hold on
    for a = 1:3
        plot(binMat{a}(:,1), binMat{a}(:,2), '.', 'color', colors(a,:), 'MarkerSize', 10);
    end
    title(sprintf('%s',beh(x).rec)); axis equal; axis square;
    xlabel('ACh minimum'); ylabel('DA maximum');
    movegui(gcf,'center');
    %}

end

% BY ANIMAL
rec = {}; for x = 1:length(beh); rec{x} = strtok(beh(x).rec,'-'); end
uni = unique(rec); nAn = length(uni);
bin_an = cell(nAn,3);
for x = 1:nAn
    ii = find(strcmp(rec,uni{x})); % index of matching recordings for this animal
    for a = 1:3
        bin_tmp = [];
        for b = 1:length(ii)
            bin_tmp = [bin_tmp; mat(ii(b)).binMat{a,1}];
        end
        bin_an{x,a} = bin_tmp;
    end
end
fprintf('DONE\n');

% NORMALIZE
bin_an_norm = cell(nAn,3); bin_dist_norm = cell(nAn,2);
mode_rew = [];
mid = edges(2:end) - bin/2;
for x = 1:nAn
    for y = 1:2
        vals_dist = [];
        for a = 1:3
            vals_raw = bin_an{x,a}(:,y);
            n = histcounts(vals_raw, edges, 'Normalization', 'probability');
            vals_dist = [vals_dist, n(:)];
        end
        [~,ii] = max(vals_dist(:,3)); mode_rew(x,y) = mid(ii); % Mode of distribution for reward
        % Normalize s.t. mode of distribution for reward is 1 and minimum
        % is 0. For ACh, mode is 0 and maximum is 1.
        bin = 0.025; edges_norm = [-3 : bin : 3]; mid_norm = edges_norm(2:end)-bin/2; % edges of normalized photometry values
        vals_dist_norm = [];
        vals_full = [bin_an{x,1}(:,y); bin_an{x,2}(:,y); bin_an{x,3}(:,y)];
        mm1 = min(vals_full); mm2 = max(vals_full);
        
        for a = 1:3
            vals_raw = bin_an{x,a}(:,y);
            switch y
                case 2
                    vals_norm = (vals_raw - mm1)./(mm2 - mm1); 
                    % vals_norm = (vals_raw - min(vals_raw))./(mode_rew(x,y) - min(vals_raw)); 
                case 1
                    vals_norm = (vals_raw - mm1)./(mm2 - mm1); 
                    % vals_norm = (vals_raw - mode_rew(x,y))./(max(vals_raw) - mode_rew(x,y)); 
            end
            vals_norm = vals_raw./abs(mode_rew(x,y)); % Normalize by mode of reward distribution
            bin_an_norm{x,a}(:,y) = vals_norm;
            n = histcounts(vals_norm, edges_norm, 'Normalization', 'probability');
            n = normalize(n, 'range');
            vals_dist_norm = [vals_dist_norm, n(:)];
        end
        bin_dist_norm{x,y} = vals_dist_norm;
    end
end

bin_dist_norm_all = cell(2,3);
for y = 1:2
    for a = 1:3
        for x = 1:nAn
            bin_dist_norm_all{y,a}(:,x) = bin_dist_norm{x,y}(:,a);
        end
    end
end

%% PLOTTING
fig = figure; fig.Position(3) = 1000;
clr = {'r','g','b'}; sm = 10;
for y = 1:2
    subplot(1,2,y); hold on
    for a = 1:3
        shadederrbar(mid_norm, movmean(nanmean(bin_dist_norm_all{y,a},2),sm), movmean(SEM(bin_dist_norm_all{y,a},2),sm), clr{a});
        % bar(mid_norm, nanmean(bin_dist_norm_all{y,a},2), 1, clr{a}, 'FaceAlpha', 0.2)
    end
    ylabel('probability'); ylim([-0.05 0.8]); axis('square');
    xlabel('photometry (norm.)');
end
subplot(1,2,1); title('ACh minima'); xlim([-3 1]);
subplot(1,2,2); title('DA maxima'); xlim([-1 3])
movegui(gcf,'center');
% normalized to mode of reward distribution maximuma/minima
% lots of overlap still!

%% STATISTICS
save_p = [];
for x = 1:nAn
    for y = 1:2
%         vec = [];
%         for a = 1:3
%             n = histcounts(bin_an{x,a}(:,y), edges, 'Normalization', 'probability');
%             vec = [vec, n(:)];
%         end
        % p = kruskalwallis(vec,[],'off');
        % [~,p] = kstest2(vec(:,1),vec(:,3));
        [~,p] = kstest2(bin_an{x,1}(:,y),bin_an{x,2}(:,y));
        % [~,p] = kstest2(bin_an_norm{x,1}(:,y),bin_an_norm{x,3}(:,y));
        save_p(x,y) = p;
    end
end

% figure; hold on
% plot([0.75 2.25],[0.05 0.05],'r');
% violinplot(save_p); xticklabels({'ACh minima','DA maxima'});
% ylabel('KS test p-value'); ylim([0 1]);

%% PLOT ANIMAL EXAMPLE
% x = 7; % used for fig2
 x = 1;% looks good

fig = figure; fig.Position([3 4]) = [1375 800]; clearvars sp
colors = [1 0 0 0.1; 0 1 0 0.1; 0 0 1 0.1];
lbl_ex = {'ACh minima','DA maxima'};
bin = 0.5;
edges = [-15 : bin : 45]; % edges of full spectrum of photometry values
mid = edges(2:end) - bin/2;
sm = 5;

for y = 1:2
    subplot(2,3,y); hold on
    % vec = [];
    for a = 1:3
        % histogram(bin_an{x,a}(:,y), 'BinEdges', edges, 'Normalization', 'probability', 'FaceColor', colors(a,[1:3]), 'FaceAlpha', 0.2, 'EdgeAlpha', 0);
        n = histcounts(bin_an{x,a}(:,y), edges, 'Normalization', 'probability');
        % vec = [vec, n(:)];
        plot(mid, movmean(n,sm), 'Color', colors(a,[1:3]));
    end
% [~,p] = kstest2(vec(:,1),vec(:,3)); [~,p(2)] = kstest2(vec(:,2),vec(:,3));
    xlabel(lbl_ex{y}); ylabel('probability'); 
    title(sprintf('%s',lbl_ex{y})); % title(sprintf('%s KS p = %1.4f, p = %1.4f',lbl_ex{y},p(1),p(2))); 
    axis square;
end
sp(4) = subplot(2,3,3); hold on
for a = 1:3
    plot(bin_an{x,a}(:,1), bin_an{x,a}(:,2), '.', 'color', colors(a,:), 'MarkerSize', 10);
    title(sprintf('%s',uni{x})); axis equal; axis square;
    xlabel(lbl_ex{1}); ylabel(lbl_ex{2});
end
for a = 1:3
    sp(a) = subplot(2,3,a+3); hold on 
    plot(bin_an{x,a}(:,1), bin_an{x,a}(:,2), '.', 'color', colors(a,:), 'MarkerSize', 10);
    xlabel(lbl_ex{1}); ylabel(lbl_ex{2}); axis square;
end; linkaxes(sp,'x'); linkaxes(sp,'y');
movegui(gcf,'center');