
%% PLOT - DA/ACh aligned to... (a) rewarded trials, (b) non-rewarded, (c) onset of lick bout
fig = figure; fig.Position(3) = 1375;
clr = {'g','m'};
for x = 1:length(lbl)
    subplot(1,3,x);
    for y = 1:2; shadederrbar(time, nanmean(align_uniAvg{x,y},2), SEM(align_uniAvg{x,y},2), clr{y}); hold on; end
    plot([0 0],[-2 16],'--k');
    xlabel('Latency (s)');
    ylabel('% dF/F');
    title(sprintf('DA/ACh to %s (n = %d)',lbl{x},length(uni))); axis square;
end
movegui(gcf,'center');

%% PLOT - Proportion of reward deliveries followed by a lick bout
a = align_uniN{1,1}; % reward yes trials for each animal
b = align_uniN{2,1}; % reward no trials
rewYesProp = a./(a+b); % proportion of trials that are rewarded per animal

figure; hold on
a = 0.8; b = 1.2;
r = a + (b-a).*rand(length(uni),1);
plot(r, rewYesProp, '.k', 'MarkerSize', 20);
errorbar(1, nanmean(rewYesProp), SEM(rewYesProp,2), '.m', 'MarkerSize', 20);
ylabel('Proportion of rewarded trials'); ylim([0 1]);
xlim([0.5 1.5]); xticks([1]); xticklabels({'Animal'}); 
title(sprintf('How many trials are rewarded (within %1.2f s)? = %1.2f +/- %1.2f',window(2), nanmean(rewYesProp), SEM(rewYesProp,2)));
axis square

%% PLOT - Timing of non-rewarded trials
rewYesWhen = [];
for x = 1:length(out)
    b = (out(x).rew_no)./length(out(x).delivery);
    rewYesWhen = [rewYesWhen; b(:)];
end

figure;
histogram(rewYesWhen,'BinWidth',0.05,'Normalization','probability');
xticks([0:0.5:1]); xticklabels({'start','mid','end'});
xlabel('non-rewarded trails w.r.t. all rewards given');
ylabel('probability');
title(sprintf('Timing of non-rewarded trials (lick within %1.2f s)?',window(2)));
 
%% When did mouse have 1st lick?
lickYes = [];
for x = 1:length(out)
    b = out(x).rew_lick - out(x).delivery(out(x).rew_yes);
    lickYes = [lickYes; b(:)];
end

figure;
histogram(lickYes,'BinWidth',0.05,'Normalization','probability');
xlabel('first lick timing w.r.t. reward (s)');
ylabel('probability');
title('When did mouse have 1st lick?');

%% Align DA and ACh to lick bouts of different vigor (duration)
    % time from first consumatory lick to last consumatory lick, which
    % should be before the next reward delivery
figure;
for x = 1:length(out)
    sp(x) = subplot(5,5,x);
lick = out(x).lick_new;
ev = out(x).delivery(out(x).rew_yes);
bin = 1/1000; window = [0 2]; % window preceding next reward delivey
peth = getClusterPETH(lick, ev, bin, window); % PETH: lick aligned to rewarded trials in 1 ms bins
cts = peth.cts{1};
[~, lick_first] = max(cts~=0, [], 1); % Find first non-zero index for each trial
[~, lick_last] = max(flipud(cts)~=0, [], 1); % Find last non-zero index
lick_last = window(2)/bin - (lick_last - 1); % Flip back because last lick was determined from flipped matrix 

% figure; 
plot(lick_last - lick_first,'.k');
hold on; 
plot([1 length(ev)], [nanmean(lick_last - lick_first) nanmean(lick_last - lick_first)], '-r','LineWidth',2);
xlabel('Reward Number'); 
ylabel('Time from 1st to Last Consumatory Lick (s)');
% title(sprintf('%s',out(x).rec));
end
linkaxes(sp,'y');
