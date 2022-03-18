%%
% fp_cell = {}; fp_max = {};
for y = 1:length(s) % iterate over infusion
    %%
    beh = s(y).s;
    [align_rew, t_rew, ev_rew] = plot_fp2event(beh,[-6 2],0);
    [align_acc, t_acc, ev_acc] = plot_fp2event(beh,[-6 2],0);
    s(y).a_rew = cell(length(beh),1);
    s(y).a_acc = cell(length(beh),1);

    a_rew_avg = []; a_acc_avg = [];
    for x = 1:length(beh) % iterate over animal
        fp = beh(x).FP{1}; Fs = 50; % extract full photometry signal
        idx_inf = [s(y).win(x,:)*60]; % infusion window
        
        ii = find(ev_rew{x} > idx_inf(1) & ev_rew{x} < idx_inf(end)); % index, wrt ev vector, where events are in infusion window
        a_rew = align_rew{x,1}(:,ii); % extract aligned traces that are in infusion window
%         a_rew = a_rew - s(y).base_mu_rest(x); % subtract mean of baseline recording to center around 0%
        a_rew = a_rew - nanmean(a_rew([find(t_rew == -1):find(t_rew == 0)],:),1);
%         m_rew = min(a_rew(find(t_rew == 0):find(t_rew == 1),:));
        
        ii = find(ev_acc{x} > idx_inf(1) & ev_acc{x} < idx_inf(end)); % index, wrt ev vector, where events are in infusion window
        a_acc = align_acc{x,1}(:,ii); % extract aligned traces that are in infusion window
%         a_acc = a_acc - s(y).base_mu_rest(x); % subtract mean of baseline recording to center around 0%
        a_acc = a_acc - nanmean(a_acc([find(t_acc == -1):find(t_acc == -0.5)],:),1);
%         m_acc = max(a_acc(find(t_acc == -0.5):find(t_acc == 0.5),:));
        
        s(y).a_rew{x} = a_rew;
        s(y).a_acc{x} = a_acc;
        a_rew_avg(:,x) = nanmean(a_rew,2); 
        a_acc_avg(:,x) = nanmean(a_acc,2);
%         fp_cell{x,y} = a_rew; fp_max{x,y} = m_rew(:);
    end
    s(y).a_rew_avg = a_rew_avg;
    s(y).a_acc_avg = a_acc_avg;
end

%% PLOT AVERAGE
fig = figure; fig.Position(3) = 1000;
clr = [0.05 0.75 0.45; 0.5 0.2 0.55; 0.85 0.35 0.1; 0.3 0.75 0.9];
subplot(1,2,1); hold on
plot([0 0],[-4 4], '--k');
for y = 1:length(s)
    shadederrbar(t_rew, nanmean(s(y).a_rew_avg,2), SEM(s(y).a_rew_avg,2), clr(y,:)); hold on
end
xlabel('Latency to Reward (s)'); xlim([-1 2]); xticks([-1:2]);
ylabel('ACh3.0 (%dF/F)'); ylim([-4 4]); yticks([-4:2:4]);
title(sprintf('ACh to reward'));
axis('square');
%
subplot(1,2,2); hold on
plot([0 0],[-2 7], '--k');
for y = 1:length(s)
    shadederrbar(t_acc, nanmean(s(y).a_acc_avg,2), SEM(s(y).a_acc_avg,2), clr(y,:)); hold on
end
xlabel('Latency to Acc Peak (s)'); xlim([-1.5 1.5]); xticks([-1:1]);
ylabel('ACh3.0 (%dF/F)'); ylim([-2 7]); yticks([-2:2:6])
title(sprintf('ACh to acceleration'));

movegui(gcf,'center');

%% Find MIN/MAX and LATENCY
nAn = 4;
amp_rew = []; lat_rew = []; amp_acc = []; lat_acc = [];
win_rew = [find(t_rew == 0):find(t_rew == 1)];
win_acc = [find(t_acc == 0):find(t_acc == 0.5)];
for y = 1:3; for x = 1:4
    [mm, ii] = min(s(y).a_rew{x}(win_rew,:)); % Pause minima within specified window
    ii = ii + win_rew(1) - 1; % Adjust latency for start time
    amp_rew(x,y) = nanmean(mm); % Average pause minima across trials
    lat_rew(x,y) = nanmean(t_rew(ii)); % Average pause latency across trials
    [mm, ii] = max(s(y).a_acc{x}(win_acc,:)); % Acc maxima within specified window
    ii = ii + win_acc(1) - 1; % Adjust latency for start time
    amp_acc(x,y) = nanmean(mm); % Average acc maxima across trials
    lat_acc(x,y) = nanmean(t_acc(ii)); 
    end; end

%% PLOT: normalize to AVERAGE aCSF
fig = figure; fig.Position([3 4]) = [1000 855];
% (1) Reward pause AMPLITUDE
a = -amp_rew./nanmean(-amp_rew(:,1)); % Normalize amplitude by aCSF
subplot(2,2,1); hold on
plot([1 3],[1 1],'--k'); % Plot line at 1 
errorbar(abs(nanmean(a)'),abs(SEM(a,1)'),'.k','MarkerSize',20); % Error bars
plot([1.15; 1.85; 2.85].*ones(3,nAn),abs(a'),'.k','MarkerSize',20); % Averages for each condition
plot([1.15; 1.85].*ones(2,nAn),abs(a(:,[1 2])'),':b','MarkerSize',20); % Connecting lines
plot([1.15; 2.85].*ones(2,nAn),abs(a(:,[1 3])'),':r','MarkerSize',20); % Connecting lines
xlim([0.5 3.5]); xticks([1:3]); xticklabels({s.inf});
ylabel('Amplitude'); ylim([0 2]); yticks([0:0.5:2]);
p = []; for x = 2:3; [~, p(x)] = ttest(a(:,1),a(:,x)); end
title(sprintf('Reward Amp (a/da: %1.3f)(a/glu: %1.3f)',p(2),p(3))); axis('square');

% (2) Reward pause LATENCY
a = lat_rew./nanmean(lat_rew(:,1)); % Normalize amplitude by aCSF
subplot(2,2,2); hold on
plot([1 3],[1 1],'--k'); % Plot line at 1 
errorbar(abs(nanmean(a)'),abs(SEM(a,1)'),'.k','MarkerSize',20); % Error bars
plot([1.15; 1.85; 2.85].*ones(3,nAn),abs(a'),'.k','MarkerSize',20); % Averages for each condition
plot([1.15; 1.85].*ones(2,nAn),abs(a(:,[1 2])'),':b','MarkerSize',20); % Connecting lines
plot([1.15; 2.85].*ones(2,nAn),abs(a(:,[1 3])'),':r','MarkerSize',20); % Connecting lines
xlim([0.5 3.5]); xticks([1:3]); xticklabels({s.inf});
ylabel('Latency'); ylim([0 3]); yticks([0:0.5:3]);
p = []; for x = 2:3; [~, p(x)] = ttest(a(:,1),a(:,x)); end
title(sprintf('Reward Lag (a/da: %1.3f)(a/glu: %1.3f)',p(2),p(3))); axis('square');

% (3) Acceleration AMPLITUDE
a = amp_acc./nanmean(amp_acc(:,1)); % Normalize amplitude by aCSF
subplot(2,2,3); hold on
plot([1 3],[1 1],'--k'); % Plot line at 1 
errorbar(abs(nanmean(a)'),abs(SEM(a,1)'),'.k','MarkerSize',20); % Error bars
plot([1.15; 1.85; 2.85].*ones(3,nAn),abs(a'),'.k','MarkerSize',20); % Averages for each condition
plot([1.15; 1.85].*ones(2,nAn),abs(a(:,[1 2])'),':b','MarkerSize',20); % Connecting lines
plot([1.15; 2.85].*ones(2,nAn),abs(a(:,[1 3])'),':r','MarkerSize',20); % Connecting lines
xlim([0.5 3.5]); xticks([1:3]); xticklabels({s.inf});
ylabel('Amplitude'); ylim([0 3]); yticks([0:0.5:3]);
p = []; for x = 2:3; [~, p(x)] = ttest(a(:,1),a(:,x)); end
title(sprintf('Acc Amp (a/da: %1.3f)(a/glu: %1.3f)',p(2),p(3))); axis('square');

% (4) Acceleration LATENCY
a = lat_acc./nanmean(lat_acc(:,1)); % Normalize amplitude by aCSF
subplot(2,2,4); hold on
plot([1 3],[1 1],'--k'); % Plot line at 1 
errorbar(abs(nanmean(a)'),abs(SEM(a,1)'),'.k','MarkerSize',20); % Error bars
plot([1.15; 1.85; 2.85].*ones(3,nAn),abs(a'),'.k','MarkerSize',20); % Averages for each condition
plot([1.15; 1.85].*ones(2,nAn),abs(a(:,[1 2])'),':b','MarkerSize',20); % Connecting lines
plot([1.15; 2.85].*ones(2,nAn),abs(a(:,[1 3])'),':r','MarkerSize',20); % Connecting lines
xlim([0.5 3.5]); xticks([1:3]); xticklabels({s.inf});
ylabel('Latency'); ylim([0 2]); yticks([0:0.5:3]);
p = []; for x = 2:3; [~, p(x)] = ttest(a(:,1),a(:,x)); end
title(sprintf('Acc Lag (a/da: %1.3f)(a/glu: %1.3f)',p(2),p(3))); axis('square');
movegui(gcf,'center');

%% PLOT: RAW VALUES
fig = figure; fig.Position([3 4]) = [1000 855];
% (1) Reward pause AMPLITUDE
a = -amp_rew; % Normalize amplitude by aCSF
subplot(2,2,1); hold on
errorbar(abs(nanmean(a)'),abs(SEM(a,1)'),'.k','MarkerSize',20); % Error bars
plot([1.15; 1.85; 2.85].*ones(3,nAn),abs(a'),'.k','MarkerSize',20); % Averages for each condition
plot([1.15; 1.85].*ones(2,nAn),abs(a(:,[1 2])'),':b','MarkerSize',20); % Connecting lines
plot([1.15; 2.85].*ones(2,nAn),abs(a(:,[1 3])'),':r','MarkerSize',20); % Connecting lines
xlim([0.5 3.5]); xticks([1:3]); xticklabels({s.inf});
ylabel('Amplitude'); % ylim([0 2]); yticks([0:0.5:2]);
p = []; for x = 2:3; [~, p(x)] = ttest(a(:,1),a(:,x)); end
title(sprintf('Reward Amp (a/da: %1.3f)(a/glu: %1.3f)',p(2),p(3))); axis('square');

% (2) Reward pause LATENCY
a = lat_rew; % Normalize amplitude by aCSF
subplot(2,2,2); hold on
errorbar(abs(nanmean(a)'),abs(SEM(a,1)'),'.k','MarkerSize',20); % Error bars
plot([1.15; 1.85; 2.85].*ones(3,nAn),abs(a'),'.k','MarkerSize',20); % Averages for each condition
plot([1.15; 1.85].*ones(2,nAn),abs(a(:,[1 2])'),':b','MarkerSize',20); % Connecting lines
plot([1.15; 2.85].*ones(2,nAn),abs(a(:,[1 3])'),':r','MarkerSize',20); % Connecting lines
xlim([0.5 3.5]); xticks([1:3]); xticklabels({s.inf});
ylabel('Latency'); % ylim([0 3]); yticks([0:0.5:3]);
p = []; for x = 2:3; [~, p(x)] = ttest(a(:,1),a(:,x)); end
title(sprintf('Reward Lag (a/da: %1.3f)(a/glu: %1.3f)',p(2),p(3))); axis('square');

% (3) Acceleration AMPLITUDE
a = amp_acc; % Normalize amplitude by aCSF
subplot(2,2,3); hold on
errorbar(abs(nanmean(a)'),abs(SEM(a,1)'),'.k','MarkerSize',20); % Error bars
plot([1.15; 1.85; 2.85].*ones(3,nAn),abs(a'),'.k','MarkerSize',20); % Averages for each condition
plot([1.15; 1.85].*ones(2,nAn),abs(a(:,[1 2])'),':b','MarkerSize',20); % Connecting lines
plot([1.15; 2.85].*ones(2,nAn),abs(a(:,[1 3])'),':r','MarkerSize',20); % Connecting lines
xlim([0.5 3.5]); xticks([1:3]); xticklabels({s.inf});
ylabel('Amplitude'); % ylim([0 3]); yticks([0:0.5:3]);
p = []; for x = 2:3; [~, p(x)] = ttest(a(:,1),a(:,x)); end
title(sprintf('Acc Amp (a/da: %1.3f)(a/glu: %1.3f)',p(2),p(3))); axis('square');

% (4) Acceleration LATENCY
a = lat_acc; % Normalize amplitude by aCSF
subplot(2,2,4); hold on
errorbar(abs(nanmean(a)'),abs(SEM(a,1)'),'.k','MarkerSize',20); % Error bars
plot([1.15; 1.85; 2.85].*ones(3,nAn),abs(a'),'.k','MarkerSize',20); % Averages for each condition
plot([1.15; 1.85].*ones(2,nAn),abs(a(:,[1 2])'),':b','MarkerSize',20); % Connecting lines
plot([1.15; 2.85].*ones(2,nAn),abs(a(:,[1 3])'),':r','MarkerSize',20); % Connecting lines
xlim([0.5 3.5]); xticks([1:3]); xticklabels({s.inf});
ylabel('Latency'); % ylim([0 2]); yticks([0:0.5:3]);
p = []; for x = 2:3; [~, p(x)] = ttest(a(:,1),a(:,x)); end
title(sprintf('Acc Lag (a/da: %1.3f)(a/glu: %1.3f)',p(2),p(3))); axis('square');
movegui(gcf,'center');
