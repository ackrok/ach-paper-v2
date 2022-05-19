%% LOAD DATA from achda_plot_CorrCoh.m
load('C:\Users\Anya\OneDrive - NYU Langone Health\TritschLab - OneDrive\Manuscripts\Anya - ACh\figures\2022Jan_figures\ach+da\achda_corr+coh_STATS_wt+b2flox.mat')
load('C:\Users\Anya\OneDrive - NYU Langone Health\TritschLab - OneDrive\Manuscripts\Anya - ACh\figures\2022Jan_figures\ach+da\achda_corr+coh_STATS_dms+dls.mat')
load('C:\Users\Anya\OneDrive - NYU Langone Health\TritschLab - OneDrive\Manuscripts\Anya - ACh\figures\2022Jan_figures\ach+da\cannula_corr+coh_STATS_acsf+dhbe.mat')

%% STATISTICS
stats = stats_dms;

x = 1; % correlation coefficient MAX
a = nanmean(stats{x}); b = SEM(stats{x},1);
fprintf('DA+ACh coeff max (R) \n   REW: %1.2f +/- %1.2f \n   MOV: %1.2f +/- %1.2f  \n   IMM: %1.2f +/- %1.2f \n',a(3),b(3),a(2),b(2),a(1),b(1))

x = 2; % correlation coefficient LAG
a = round(nanmean(stats{x})); b = round(SEM(stats{x},1));
fprintf('DA+ACh lag (ms) \n   REW: %d +/- %d ms \n   MOV: %d +/- %d ms \n   IMM: %d +/- %d ms \n',a(3),b(3),a(2),b(2),a(1),b(1))

x = 3; % coherence magnitude
a = nanmean(stats{x}); b = SEM(stats{x},1);
fprintf('DA+ACh coherence magnitude \n   REW: %1.2f +/- %1.2f \n   MOV: %1.2f +/- %1.2f  \n   IMM: %1.2f +/- %1.2f \n',a(3),b(3),a(2),b(2),a(1),b(1))

x = 4; % phase offset
a = round(nanmean(stats{x})); b = round(SEM(stats{x},1));
fprintf('DA+ACh phase offset (degrees) \n   REW: %d +/- %d degrees\n   MOV: %d +/- %d degrees \n   IMM: %d +/- %d degrees \n n = %d mice \n',a(3),b(3),a(2),b(2),a(1),b(1),size(stats{x},1))

p = []; for x = 1:4; p(x) = kruskalwallis(stats{x},[],'off'); end
fprintf('p-value for KS test comparing IMM, MOV, REW \n   coeff: %1.4f , coeff lag: %1.4f \n   coh mag: %1.4f , phase: %1.4f \n',p(1),p(2),p(3),p(4))

p = []; for x = 1:4; p(x) = anova1(stats{x},[],'off'); end
fprintf('p-value for ANOVA comparing IMM, MOV, REW \n   coeff: %1.4f , coeff lag: %1.4f \n   coh mag: %1.4f , phase: %1.4f \n',p(1),p(2),p(3),p(4))

%% PLOT/STATS: WT vs b2flox
stats_1 = stats_wt; stats_2 = stats_b2; % CHANGE

fig = figure; fig.Position([3 4]) = [1000 850]; hold on; movegui(gcf,'center');
for x = 1:4
    subplot(2,2,x); hold on;
    a = [stats_1{x}, [stats_2{x};nan(7,3)]]; a = a(:,[1 4 2 5 3 6]);
    violinplot(a); % violinplot values
    errorbar(0.25+[1:size(a,2)], nanmean(a), SEM(a,1), '.k', 'MarkerSize', 20); % errorbar
    xticklabels({'imm','imm-b2','mov','mov-b2','rew','rew-b2'});
    ylabel(sprintf('%s',lbl{x})); ylim([lims(x,:)]);
    p = []; for b = 1:3; p(b) = ranksum(stats_1{x}(:,b), stats_2{x}(:,b)); end
    title(sprintf('ranksum imm=%1.3f | mov=%1.3f | rew=%1.3f',p(1),p(2),p(3))); axis('square');
end
% x = 1; subplot(2,2,x); ylim([-1 0]); yticks([-1:0.2:0]); % adjust axes
% x = 2; subplot(2,2,x); ylim([-250 0]); yticks([-250:50:0]); % adjust axes
% x = 3; subplot(2,2,x); ylim([0 1]); yticks([0:0.2:1]); % adjust axes
x = 4; subplot(2,2,x); ylim([-180 0]); yticks([-180:45:180]); % adjust axes

%% PLOT/STATS: aCSF vs nAChR antag
stats_1 = stats_acsf; stats_2 = stats_dhbe; % CHANGE

fig = figure; fig.Position([3 4]) = [1000 850]; hold on; movegui(gcf,'center');
for x = 1:4
    subplot(2,2,x); hold on;
    a = [stats_1{x}, stats_2{x}]; a = a(:,[1 4 2 5 3 6]); a(3,:) = nan;
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

%% PLOT/STATS: WT vs b2flox (v2, % control)
stats_1 = stats_wt; stats_2 = stats_b2; % CHANGE

fig = figure; fig.Position([3 4]) = [1000 850]; hold on; movegui(gcf,'center');
for x = 1:4
    subplot(2,2,x); hold on;
    a = stats_2{x}./nanmean(stats_1{x});
    plot([0.5 3.5],[1 1],'-k');
    violinplot(a); % violinplot values
%     plot(a','--','Color',[0 0 0 0.2]) % connect values
    errorbar([1.25 2.25 3.25], nanmean(a), SEM(a,1), '.-k', 'MarkerSize', 20); % errorbar
    xticklabels({'imm','mov','rew'});
    ylabel(sprintf('%s',lbl{x})); ylim([0 1.5]); yticks([0:0.25:1.5]);
    % p = kruskalwallis(a,[],'off'); title(sprintf('kruskalwallis p=%1.3f',p));
    % p = anova1(a,[],'off'); title(sprintf('anova p=%1.3f',p));
    p2 = []; for b = 1:3; [~,p2(b)] = ttest2(stats_1{x}(:,b), stats_2{x}(:,b)); end
    title(sprintf('ttest imm=%1.3f | mov=%1.3f | rew=%1.3f',p2(1),p2(2),p2(3))); axis('square');
    axis('square');
end

%% STATS: DMS vs DLS
x = 1; [~,~,stats] = anova1(coh_dmsdls{x},[],'off');
c = multcompare(stats,'display','off');
fprintf('p-value for ANOVA - coherence - immobility \n   DA-DLS/ACh-DMS: %1.4f , DA-DMS/ACh-DLS: %1.4f vs DA-DLS/ACh-DLS \n',c(2,6),c(3,6))

x = 2; [~,~,stats] = anova1(coh_dmsdls{x},[],'off');
c = multcompare(stats,'display','off');
fprintf('p-value for ANOVA - coherence - locomotion \n   DA-DLS/ACh-DMS: %1.4f , DA-DMS/ACh-DLS: %1.4f vs DA-DLS/ACh-DLS \n',c(2,6),c(3,6))

x = 3; [~,~,stats] = anova1(coh_dmsdls{x},[],'off');
c = multcompare(stats,'display','off');
fprintf('p-value for ANOVA - coherence - reward \n   DA-DLS/ACh-DMS: %1.4f , DA-DMS/ACh-DLS: %1.4f vs DA-DLS/ACh-DLS \n',c(2,6),c(3,6))