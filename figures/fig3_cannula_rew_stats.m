%%
load('C:\Users\Anya\Desktop\FP_LOCAL\AK189-197_cannula+rew_v4.mat')
% fp_cell = {}; fp_max = {};
for y = 1:length(s) % iterate over infusion
    %%
    beh = s(y).s;
    [align_rew, t_rew, ev_rew] = plot_fp2event(beh,[-6 2],0);
%     [align_acc, t_acc, ev_acc] = plot_fp2event(beh,[-6 2],0);
    s(y).a_rew = cell(length(beh),1);
%     s(y).a_acc = cell(length(beh),1);

    a_rew_avg = []; % a_acc_avg = [];
    for x = 1:length(beh) % iterate over animal
        fp = beh(x).FP{1}; Fs = 50; % extract full photometry signal
        idx_inf = [s(y).win(x,:)*60]; % infusion window
        
        ii = find(ev_rew{x} > idx_inf(1) & ev_rew{x} < idx_inf(end)); % index, wrt ev vector, where events are in infusion window
        a_rew = align_rew{x,1}(:,ii); % extract aligned traces that are in infusion window
%         a_rew = a_rew - s(y).base_mu_rest(x); % subtract mean of baseline recording to center around 0%
        a_rew = a_rew - nanmean(a_rew([find(t_rew == -1):find(t_rew == 0)],:),1);
%         m_rew = min(a_rew(find(t_rew == 0):find(t_rew == 1),:));
        
%         ii = find(ev_acc{x} > idx_inf(1) & ev_acc{x} < idx_inf(end)); % index, wrt ev vector, where events are in infusion window
%         a_acc = align_acc{x,1}(:,ii); % extract aligned traces that are in infusion window
% %         a_acc = a_acc - s(y).base_mu_rest(x); % subtract mean of baseline recording to center around 0%
%         a_acc = a_acc - nanmean(a_acc([find(t_acc == -1):find(t_acc == -0.5)],:),1);
% %         m_acc = max(a_acc(find(t_acc == -0.5):find(t_acc == 0.5),:));
        
        s(y).a_rew{x} = a_rew;
%         s(y).a_acc{x} = a_acc;
        a_rew_avg(:,x) = nanmean(a_rew,2); 
%         a_acc_avg(:,x) = nanmean(a_acc,2);
%         fp_cell{x,y} = a_rew; fp_max{x,y} = m_rew(:);
    end
    s(y).a_rew_avg = a_rew_avg;
%     s(y).a_acc_avg = a_acc_avg;
end

%% PLOT EXAMPLE
fig = figure; %fig.Position(3) = 1000;
clr = [0.05 0.75 0.45; 0.5 0.2 0.55; 0.85 0.35 0.1; 0.3 0.75 0.9];

x = 3;
%subplot(1,2,1); hold on
plot([0 0],[-4 4], '--k');
for y = 1:3
    shadederrbar(t_rew, nanmean(s(y).a_rew{x},2), SEM(s(y).a_rew{x},2), clr(y,:)); hold on
end
xlabel('Latency to Reward (s)'); xlim([-1 2]); xticks([-1:2]);
ylabel('ACh3.0 (%dF/F)'); ylim([-4 4]); yticks([-4:2:4]);
title(sprintf('ACh to reward'));
axis('square');

%%
% Find MIN/MAX and LATENCY
nAn = 4;
amp_rew = []; lat_rew = [];
win_rew = [find(t_rew == 0):find(t_rew == 1)];
win_base = [find(t_rew == -6):find(t_rew == -1)];

win_new = [];
amp_adj = []; amp_base = []; amp_base_low = [];
amp_new = []; amp_low = []; amp_upp = [];
lat_new = [];

for y = 1 for x = 1:4
    [mm, ii] = min(s(y).a_rew{x}(win_rew,:)); % Pause minima within specified window
    ii = ii + win_rew(1) - 1; % Adjust latency for start time
    prc = prctile(ii, [40 50 60]);
    win_new(x,:) = round(prc);
    amp_rew(x,y) = nanmean(mm); % Average pause minima across trials
    
    end; end

for y = 1:3 for x = 1:4
    mm = s(y).a_rew{x}([win_new(x,1):win_new(x,3)],:); % Pause minima within specified window
    mm(mm > 0) = nan;
    amp_new(x,y) = nanmean(nanmean(mm));
    amp_low(x,y) = nanmean(SEM(mm,1)); amp_upp(x,y) = nanmean(SEM(mm,1));
    % amp_new(x,y) = nanmean(prctile(mm,50)); 
    % amp_low(x,y) = nanmean(prctile(mm,50))-nanmean(prctile(mm,25)); amp_upp(x,y) = nanmean(prctile(mm,75))-nanmean(prctile(mm,50));
    
    win_base = win_new(x,:) - 2.*(win_new(x,:) - find(t_rew == 0));
    mm_base = s(y).a_rew{x}([win_base(3):win_base(1)],:); % Pause minima within specified window
    amp_base(x,y) = nanmean(nanmean(mm_base));
    amp_base_low(x,y) = nanmean(SEM(mm_base,1));
    % amp_base(x,y) = nanmean(prctile(mm_base,50));
    % amp_base_low(x,y) = nanmean(prctile(mm_base,50))-nanmean(prctile(mm_base,25)); 
    % amp_base_upp(x,y) = nanmean(prctile(mm_base,75))-nanmean(prctile(mm_base,50));
    
    amp_adj(x,y) = amp_new(x,y) - amp_base(x,y);
    end; end  

%% PLOT FINAL
figure; hold on
for y = 1:3
    errorbar(y.*ones(nAn,1), amp_adj(:,y), amp_low(:,y), amp_upp(:,y), '.', 'MarkerSize',20,'Color',clr(y,:));
end
plot([1;2].*ones(2,nAn),amp_adj(:,[1 2])','-','Color',[0 0 0 0.2]);
plot([1;3].*ones(2,nAn),amp_adj(:,[1 3])','-','Color',[0 0 0 0.2]);
errorbar([1:3]-0.25, nanmean(amp_adj,1), SEM(amp_adj,1), '.k', 'MarkerSize',20);
xticks([1:3]); xlim([0.5 3.5]);xticklabels({'acsf','d1d2','glu'});
ylabel('avg %dF/F within pause window - baseline'); ylim([-5 0]); yticks([-5:1]);
% [p,~,stats] = kruskalwallis(amp_adj,[],'off');
[p,~,stats] = anova1(amp_adj,[],'off');
c = multcompare(stats,'display','off');
% title(sprintf('ks=%1.3f \n acsf/d1d2=%1.3f | acsf/glu=%1.3f',p,c(1,6),c(2,6))); axis('square');
title(sprintf('anova=%1.3f \n acsf/d1d2=%1.3f | acsf/glu=%1.3f',p,c(1,6),c(2,6))); axis('square');

%% PLOT + BASELINE
figure; hold on
errorbar(t_rew(win_new(:,2)).*1000,amp_base(:,1), amp_base(:,2), '.k', 'MarkerSize',20);
for y = 1:3
    errorbar(t_rew(win_new(:,2)).*1000, amp_new(:,y), amp_low(:,y), amp_upp(:,y), '.', 'MarkerSize',20,'Color',clr(y,:));
end
legend({'base','acsf','d1d2','glu'});
xlabel('time to reward (ms)'); xlim([0 600]); xticks([0:100:1000]); 
ylabel('avg %dF/F within pause window'); ylim([-5 1]); yticks([-5:1]);
title(sprintf('ks=0.0051 \n base/acsf=0.0098 | base/d1d2=0.0247 | rest n.s.')); axis('square');

 %% Find MIN/MAX and LATENCY
% nAn = 4;
% amp_rew = []; lat_rew = [];
% win_rew = [find(t_rew == 0):find(t_rew == 1)];
% win_base = [find(t_rew == -6):find(t_rew == -1)];
% 
% win_new = [];
% win_base = []; amp_base = [];
% amp_new = []; amp_low = []; amp_upp = [];
% 
% for y = 1 for x = 1:4
%     [mm, ii] = min(s(y).a_rew{x}(win_rew,:)); % Pause minima within specified window
%     ii = ii + win_rew(1) - 1; % Adjust latency for start time
    % base = s(y).a_rew{x}(win_base,:);
    % base = 1.5*nanmean(nanstd(base)); % numStd = 2
    % if y == 3; mm(abs(mm) < base) = 0; ii(abs(mm) < base) = 50; end

%   subplot(1,2,2); hold on  
%     ii = ii + win_rew(1) - 1; % Adjust latency for start time
%     ii_ms = t_rew(ii)*1000;
%     
%     errorbar(nanmean(ii_ms),nanmean(-mm),SEM(mm,2),SEM(mm,2),SEM(ii_ms,2),SEM(ii_ms,2),'.','MarkerSize',20,'Color',clr(y,:));
%     amp_rew(x,y) = nanmean(mm); % Average pause minima across trials
%     lat_rew(x,y) = nanmean(t_rew(ii)); % Average pause latency across trials
%     end; end
% xlabel('time from reward (ms)'); xlim([0 1000]);
% ylabel('reward pause amplitude'); ylim([0 5]); yticks([0:5]);
% title('reward pause'); axis('square');

%% LAG/VAL at PAUSE MAXIMUM for aCSF/D1D2 comparison
lag = []; val = [];
r = [0.2 0.7]; r = (r/(1/Fs) + find(t_rew==0)); % range for ACh trough
for y = 1:2; for x = 1:4
    [a,b] = min(s(y).a_rew{x}(r(1):r(2),:)); % find local minima within range
    c = t_rew(b + r(1) - 1); c(isnan(a)) = nan; % convert index to seconds
    lag(x,y) = nanmean(c);
    val(x,y) = nanmean(a);
    end; end

figure; hold on
y = 1; plot(lag(:,y).*1000, val(:,y), '.k', 'MarkerSize', 20);
y = 2; plot(lag(:,y).*1000, val(:,y), '.g', 'MarkerSize', 20);
y = 1; errorbar(nanmean(lag(:,y).*1000), nanmean(val(:,y)), SEM(val(:,y),1), SEM(val(:,y),1), SEM(lag(:,y).*1000,1), SEM(lag(:,y).*1000,1), '.k', 'MarkerSize', 20);
y = 2; errorbar(nanmean(lag(:,y).*1000), nanmean(val(:,y)), SEM(val(:,y),1), SEM(val(:,y),1), SEM(lag(:,y).*1000,1), SEM(lag(:,y).*1000,1), '.g', 'MarkerSize', 20);
ylabel('ACh trough amp'); ylim([-5 0]); yticks([-8:2:0]);
xlabel('time to rew (s)'); xlim([200 700]);
[~,p] = ttest2(lag(:,1),lag(:,2)); [~,p(2)] = ttest2(val(:,1),val(:,2));
title(sprintf('ttest2: lag = %1.3f, val = %1.3f',p(1),p(2))); axis('square');