%% 
load('C:\Users\Anya\Desktop\FP_LOCAL\AK189-197_cannula+rew_v4.mat')
cannula_immPausePeak
cannula_RewOnsetPausePeak

%%
ach2ach_acsf = cell(2,2); 
dur = []; amp = [];
dur_imm = []; amp_imm = []; freq_imm = [];

y = 1; % aCSF infusion
for x = 1:4 % iteratre over animals
    for z = 1:2 % iterate over pause/peak
        a = s(y).ach2ach{x,z};
        a = a - nanmean(a(find(sta_time == -4):find(sta_time == -1),:)); % subtract baseline [-4 -1]
        a = a(find(sta_time == -1):find(sta_time == 1),:); % restrict to [-1 1] window
        a = nanmean(a,2);
        switch z; case 1; b = (a - min(a))/(a(1,:) - min(a)); % min max normalization
            case 2; b = (a - a(1,:))/(max(a) - a(1,:)); end
        ach2ach_acsf{z,1}(:,x) = b;
        
        freq_imm(x,z) = s(y).freq(x,z);
        dur_imm(x,z) = nanmean(s(y).dur{x,z}); 
        amp_imm(x,z) = nanmean(s(y).amp{x,z});
        
        %%
        if z == 2; if x == 2
                ach2ach_acsf{z,2}(:,x) = nan(101,1); dur(x,z) = 12; amp(x,z) = 5.1255;
                continue; end; end
        a2 = s(y).ach2rew{x,z};
        a2 = a2 - nanmean(a2(find(sta_time == -4):find(sta_time == -1),:)); % subtract baseline [-4 -1]
        a2 = a2(find(sta_time == -1):find(sta_time == 1),:); % restrict to [-1 1] window
        a2 = nanmean(a2,2);
        switch z; case 1; b2 = (a2 - min(a2))/(a2(1,:) - min(a2)); % min max normalization
            case 2; b2 = (a2 - a2(1,:))/(max(a2) - a2(1,:)); end
        ach2ach_acsf{z,2}(:,x) = b2;
        
        sig = a2;
        maxVal = sig(51);
        switch z
            case 1
                tmp = find(sig < 0.5*maxVal); % Find index at half max
                halfMax = tmp(1); % First value that crosses threshold 
                tmp = find(sig(51:end) > 0.5*maxVal); % Find locations that exceed new threshold, after maximum deflection point
                halfMax(2) = tmp(1)+50;
            case 2
                tmp = find(sig > 0.5*maxVal); % Find index at half max
                halfMax = tmp(1); % First value that crosses threshold 
                tmp = find(sig(51:end) < 0.5*maxVal); % Find locations that exceed new threshold, after maximum deflection point
                halfMax(2) = tmp(1)+50;
        end
        dur(x,z) = diff(halfMax);
        amp(x,z) = maxVal;
    end
end
t_sub = sta_time(find(sta_time == -1):find(sta_time == 1));

%%
fig = figure; fig.Position(3) = 1000;
for z = 1:2
    switch z; case 1; lbl = 'pause'; case 2; lbl = 'peak'; end
    subplot(1,2,z);
    % plot(t_sub, ach2ach_acsf{z,1});
    plot(t_sub, ach2ach_acsf{z,2}); legend
    % plot(t_sub, ach2ach_acsf{z,1}, 'Color', [0.05 0.75 0.45 0.2]);
    % plot(t_sub, ach2ach_acsf{z,2}, 'Color', [0 0 0 0.2]);
    % shadederrbar(t_sub, nanmean(ach2ach_acsf{z,1},2), SEM(ach2ach_acsf{z,1},2), [0.05 0.75 0.45]);
    % shadederrbar(t_sub, nanmean(ach2ach_acsf{z,2},2), SEM(ach2ach_acsf{z,2},2), 'k');
    xlabel('latency (s)'); xlim([-1 1]);
    ylabel('ACh (normalized)');
    title(sprintf('ACh to ACh %s - aCSF',lbl)); axis('square');
end

%%
z = 2;
switch z; case 1; lbl = 'Pause'; case 2; lbl = 'Peak'; end

fig = figure; fig.Position(3) = 1375;
subplot(1,3,1); hold on
shadederrbar(t_sub, nanmean(ach2ach_acsf{z,1},2), SEM(ach2ach_acsf{z,1},2), [0.05 0.75 0.45]);
shadederrbar(t_sub, nanmean(ach2ach_acsf{z,2},2), SEM(ach2ach_acsf{z,2},2), 'k');
xlabel('latency (s)'); xlim([-1 1]);
ylabel('ACh (normalized)');
title(sprintf('%s. freq = %1.2f +/- %1.2f Hz',lbl,nanmean(freq_imm(:,z)),SEM(freq_imm(:,z),1))); axis('square');

subplot(1,3,2); hold on
plot_me = [amp(:,z), amp_imm(:,z)]; % CHANGE
errorbar(abs(nanmean(plot_me)'),abs(SEM(plot_me,1)'),'.k','MarkerSize',20);
plot([1.15; 1.85].*ones(2,4),abs([plot_me']),'.:b','MarkerSize',20);
ylabel('amplitude (%dF/F)'); 
switch z; case 1; ylim([2 7]); yticks([2:7]); case 2; ylim([0 15]); yticks([0:5:15]); end
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'rew','imm'});
[~,p] = ttest(plot_me(:,1),plot_me(:,2));
p(2) = signrank(plot_me(:,1),plot_me(:,2));
title(sprintf('%s (ttest = %1.4f) (sr = %1.4f)',lbl,p(1),p(2))); axis('square')

subplot(1,3,3); hold on
plot_me = [dur(:,z), dur_imm(:,z)].*(1000/50); % CHANGE
errorbar(nanmean(plot_me)',SEM(plot_me,1)','.k','MarkerSize',20);
plot([1.15; 1.85].*ones(2,4),[plot_me'],'.:b','MarkerSize',20);
ylabel('duration (ms)');
switch z; case 1; ylim([0 400]); yticks([0:100:500]); case 2; ylim([0 500]); yticks([0:100:500]); end
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'rew','imm'});
[~,p] = ttest(plot_me(:,1),plot_me(:,2));
p(2) = signrank(plot_me(:,1),plot_me(:,2));
title(sprintf('%s (ttest = %1.4f) (sr = %1.4f)',lbl,p(1),p(2))); axis('square')

movegui(fig,'center');