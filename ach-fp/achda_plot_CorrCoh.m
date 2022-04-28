%% bzCoherogram
% need: coher_an, phase_an

%% achda_corr
% need: min_val, min_lag

%% PLOTTING
r = [6:42]; [~,r2] = min(abs(f - 2)); % range for 0.5-4Hz
coher_avg = []; phase_avg = [];
for y = 1:length(coher_an)
    % coher_avg(:,y) = coher_an{y}(r2,:); phase_avg(:,y) = phase_an{y}(r2,:); % value at 2Hz
    % coher_avg(:,y) = nanmean(coher_an{y}(r,:)); phase_avg(:,y) = nanmean(phase_an{y}(r,:)); % average within frequency band
    coher_avg(:,y) = median(coher_an{y}(r,:)); phase_avg(:,y) = median(phase_an{y}(r,:)); % median within frequency band
end

stats = cell(4,1); % initialize
lbl = {'coeff-mag','coeff-lag (ms)','coh-mag','coh-ph (deg)'};
lims = [-1 0; -250 0; 0 1; -180 0];
stats{1} = min_val; stats{2} = min_lag.*1000;
stats{3} = coher_avg; stats{4} = rad2deg(phase_avg);

%%
fig = figure; fig.Position([3 4]) = [1000 850]; hold on; movegui(gcf,'center');
for x = 1:4
    subplot(2,2,x); hold on;
    a = stats{x};
    violinplot(a); % violinplot values
%     plot(a','--','Color',[0 0 0 0.2]) % connect values
    errorbar([1.25 2.25 3.25], nanmean(a), SEM(a,1), '.-k', 'MarkerSize', 20); % errorbar
    xticklabels({'imm','mov','rew'});
    ylabel(sprintf('%s',lbl{x})); ylim([lims(x,:)]);
    p = kruskalwallis(a,[],'off');
    title(sprintf('kruskalwallis p=%1.3f',p)); axis('square');
end

% x = 1; subplot(2,2,x); ylim([-1 0]); yticks([-1:0.2:0]); % adjust axes
% x = 2; subplot(2,2,x); ylim([-250 0]); yticks([-250:50:0]); % adjust axes
% x = 3; subplot(2,2,x); ylim([0 1]); yticks([0:0.2:1]); % adjust axes
x = 4; subplot(2,2,x); ylim([-180 0]); yticks([-180:45:180]); % adjust axes

%% control VS b2flox stats
p = [];
for a = 1:4; for b = 1:3
     p(a,b) = ranksum(stats_wt{a}(:,b),stats_b2{a}(:,b));
%    p(a,b) = signrank(stats_acsf{a}(:,b), stats_dhbe{a}(:,b)); 
    % [~,p(a,b)] = ttest(stats_acsf{a}(:,b), stats_dhbe{a}(:,b)); 
    end; end

%% Kruskal-Wallis test
stats = stats_wt;
p = []; for a = 1:4; p(a) = kruskalwallis(stats{a},[],'off'); end

% [p, ~, stats] = kruskalwallis(stats{a})
% multcompare(stats)

%% aCSF vs DHbE
stats_1 = stats_acsf; stats_2 = stats_dhbe; % CHANGE

fig = figure; fig.Position([3 4]) = [1000 850]; hold on; movegui(gcf,'center');
for x = 1:4
    subplot(2,2,x); hold on;
    a = [stats_1{x}, stats_2{x}]; a = a(:,[1 4 2 5 3 6]);
    violinplot(a); % violinplot values
    plot(a','--','Color',[0 0 0 0.2]) % connect values
    errorbar(0.25+[1:size(a,2)], nanmean(a), SEM(a,1), '.k', 'MarkerSize', 20); % errorbar
    xticklabels({'imm-A','imm-N','mov-A','mov-N','rew-A','rew-N'});
    ylabel(sprintf('%s',lbl{x})); ylim([lims(x,:)]);
    % p = kruskalwallis(a,[],'off');
    p2 = []; for b = 1:3; [~,p2(b)] = ttest(stats_1{x}(:,b), stats_2{x}(:,b)); end
    title(sprintf('ttest imm=%1.3f | mov=%1.3f | rew=%1.3f',p2(1),p2(2),p2(3))); axis('square');
end

% x = 1; subplot(2,2,x); ylim([-1 0]); yticks([-1:0.2:0]); % adjust axes
% x = 2; subplot(2,2,x); ylim([-250 0]); yticks([-250:50:0]); % adjust axes
% x = 3; subplot(2,2,x); ylim([0 1]); yticks([0:0.2:1]); % adjust axes
x = 4; subplot(2,2,x); ylim([-180 0]); yticks([-180:45:180]); % adjust axes
