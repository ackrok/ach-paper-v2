acc2acc = cell(4,4); acc_avg = cell(4,1);
lick2rew = cell(4,4); lick_avg = cell(4,1); lick_rate = [];
for y = 1:length(s)
beh = s(y).s;
    for x = 1:length(beh)
        %% acceleration
        Fs = beh(x).Fs;
        acc = getAcc(beh(x).vel); % Extract acceleration signal
        [~,locs] = findpeaks(acc,'MinPeakProminence',0.5,'MinPeakDistance',0.5); % Location of peaks, using findpeaks function
        ev = beh(x).time(locs); % Convert peak locations to seconds
        [sta,time] = getSTA(acc, ev, Fs, [-1 1]); % Align acceleration to acceleration locs
        acc2acc{x,y} = sta;
        acc_avg{y}(:,x) = nanmean(sta,2);
        %% reward
        if all(logical(~rem(beh(x).on,1))); diffFs = 1; else; diffFs = 50; end
        if isempty(beh(x).reward); continue; end
        rew = beh(x).reward/(Fs/diffFs); % Extract reward delivery times, adjusting event times to be in seconds
        lick = beh(x).lick/(Fs/diffFs);
        bin = 1/1000; window = [0 1];
        peth = getClusterPETH(lick, rew, bin, window); % PETH: lick aligned to reward in 1 ms bins
        cts = peth.cts{1}; % Lick counts in 1ms bins for each reward trial
        lick_rate(x,y) = nanmean(sum(cts,1));
        rew_lick0 = find(sum(cts,1) == 0); % Find reward index where total licks within window is 0
        bin = 1/1000; window = [0 0];
        peth = getClusterPETH(lick, rew, bin, window);
        rew_prewlick = find(sum(peth.cts{1},1) >= 1); % Find reward index for trials where mouse licks preceding reward
        rew([rew_lick0, rew_prewlick]) = nan; 
        cts(:, [rew_lick0, rew_prewlick]) = nan; % Remove non-rewarded trials and trials where mouse licks preceding reward
        ev = rew; % Event time is rewarded trials
        [sta,time] = getSTA(beh(x).lickVec, ev, Fs, [-1 1]); % Align acceleration to acceleration locs
        lick2rew{x,y} = sta;
        lick_avg{y}(:,x) = nanmean(sta,2);
    end
end
fprintf('Done! \n');
%% average across animals
fig = figure; fig.Position(3) = 1000;
subplot(1,2,1); hold on
plot([0 0],[-0.3 0.6],'--k');
y = 1; shadederrbar(time, nanmean(acc_avg{y},2), SEM(acc_avg{y},2), 'k');
y = 3; shadederrbar(time, nanmean(acc_avg{y},2), SEM(acc_avg{y},2), 'r');
xlabel('acc peak (s)'); ylabel('acceleration');
title('acc n = 4'); axis('square');

subplot(1,2,2); hold on
plot([0 0],[0 0.6],'--k');
y = 1; shadederrbar(time, nanmean(lick_avg{y},2), SEM(lick_avg{y},2), 'k');
y = 3; shadederrbar(time, nanmean(lick_avg{y},2), SEM(lick_avg{y},2), 'r');
xlabel('reward (s)'); ylabel('lick');
title('lick n = 4'); axis('square');

%% by animal
figure;
for x = 1:4
sp(x) = subplot(2,2,x); hold on
y = 1; shadederrbar(time, nanmean(acc2acc{x,y},2), SEM(acc2acc{x,y},2), 'k');
y = 3; shadederrbar(time, nanmean(acc2acc{x,y},2), SEM(acc2acc{x,y},2), 'r');
xlabel('acc peak (s)'); ylabel('acceleration');
title(sprintf('%s',strtok(s(1).s(x).rec,'_')));
end
linkaxes(sp,'y');

figure;
for x = 1:4
sp(x) = subplot(2,2,x); hold on
y = 1; shadederrbar(time, nanmean(lick2rew{x,y},2), SEM(lick2rew{x,y},2), 'k');
y = 3; shadederrbar(time, nanmean(lick2rew{x,y},2), SEM(lick2rew{x,y},2), 'r');
xlabel('reward (s)'); ylabel('lick');
title(sprintf('%s',strtok(s(1).s(x).rec,'_')));
end
linkaxes(sp,'y');
