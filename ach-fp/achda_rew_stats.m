
%% PLOT EXAMPLE
fig = figure; %fig.Position(3) = 1000;
clr = [0.05 0.75 0.45; 0.5 0.2 0.55; 0.85 0.35 0.1; 0.3 0.75 0.9];

x =1;
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
xticks([1:3]); xlim([0.5 3.5]);xticklabels({'acsf','d1d2','glu'});
ylabel('avg %dF/F within pause window - baseline'); ylim([-5 0]); yticks([-5:1]);
[p,~,stats] = kruskalwallis(amp_adj,[],'off');
c = multcompare(stats,'display','off');
title(sprintf('ks=%1.3f \n acsf/d1d2=%1.3f | acsf/glu=%1.3f',p,c(1,6),c(2,6))); axis('square');

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