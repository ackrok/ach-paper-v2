%% RUN cannula_immPausePeak
%
%%
da2ach = cell(4,2); ach2ach = cell(4,2);
for x = 1:4
    da2ach{x,1} = s(1).da2ach{x+4,2}; % Extract DA aligned to ACh peaks for this infusion condition
    ach2ach{x,1} = s(1).ach2ach{x+4,2}; % Extract ACh aligned to ACh peaks
    if x == 3
        da2ach{x,2} = nan(401,1); % Correct for no recording with AK193
        ach2ach{x,2} = nan(401,1); 
    elseif x == 4
        da2ach{x,2} = s(4).da2ach{3,2}; ach2ach{x,2} = s(4).ach2ach{3,2};
    else
        da2ach{x,2} = s(4).da2ach{x,2}; ach2ach{x,2} = s(4).ach2ach{x,2}; % Extract DA aligned to ACh peaks for this infusion condition
    end
end
da2ach_avg = cell(2,1);
da2ach_full = cell(2,4); ach2ach_full = cell(2,1);
da2ach_norm = cell(2,1);
time_seg = sta_time(find(sta_time == -1):find(sta_time == 1));
for y = 1:2; for x = 1:4
        da_base = da2ach{x,y}(find(sta_time == -6):find(sta_time == -1),:);
        ach_base = ach2ach{x,y}(find(sta_time == -6):find(sta_time == -1),:);
        
        da_seg = da2ach{x,y}(find(sta_time == -1):find(sta_time == 1),:);
        da_seg = da_seg - nanmean(da_base,1);% da_seg(26,:); %  % Adjust for baseline
        ach_seg = ach2ach{x,y}(find(sta_time == -1):find(sta_time == 1),:);
        ach_seg = ach_seg - nanmean(ach_base,1); % Adjust for baseline
        
        da2ach_avg{y}(:,x) = nanmean(da_seg,2); % Extract average across all ACh peaks
        da2ach_full{y,x} = da_seg; % Extract DA aligned to all ACh peaks
        ach2ach_full{y} = [ach2ach_full{y}, ach_seg];
        da2ach_norm{y}(:,x) = normalize(da2ach_avg{y}(:,x),'range');
    end; end


% t_acsf = [0.16 0.18 0.02 0.04]; % t_acsf = [0.16 0.18 0.14 0.22]; 
% v_acsf = da2ach_avg{1}(time_seg == t_acsf);
% t_dhbe = [0.08 0.14 nan 0.1]; %t_dhbe = [0.08 0.14 nan 0.1]; 
% v_dhbe = da2ach_avg{2}(time_seg == t_dhbe); v_dhbe = [v_dhbe(1:2);nan;v_dhbe(3)];
% t_delta = [t_acsf(:), t_dhbe(:)]; t_delta = t_delta([1 2 4],:); 
% t_delta = t_delta./nanmean(t_delta(:,1));
% v_delta = [v_acsf(:), v_dhbe(:)]; v_delta = v_delta([1 2 4],:); 
% v_delta = v_delta./nanmean(v_delta(:,1));

v_diff = []; t_diff = [];
v_max = []; t_max = [];
for y = 1:length(da2ach_avg); for x = 1:size(da2ach_avg{y},2)
        [a, ii1] = min(da2ach_avg{y}(find(time_seg == -0.26):find(time_seg == 0),x)); % find MIN within range preceding ACh peak
        [b, ii2] = max(da2ach_avg{y}(find(time_seg == 0):find(time_seg == 0.16),x)); % find MAX within range following ACh peak
        v_diff(x,y) = b - a;
        t_diff(x,y) = time_seg(ii2 + find(time_seg == 0) - 1) - time_seg(ii1 + find(time_seg == -0.26) - 1);
        v_max(x,y) = b; t_max(x,y) = time_seg(ii2 + find(time_seg == 0) - 1);
    end
end

%% ACTUAL VALUES
fig = figure; fig.Position([3 4]) = [1000 900];
nAn = 3;
subplot(2,2,1); hold on
a = t_max([1 2 4],:).*1000;
errorbar(nanmean(a)',SEM(a,1)','.k','MarkerSize',20);
plot([1.15; 1.85].*ones(2,nAn),a','.k','MarkerSize',20);
plot([1.15; 1.85].*ones(2,nAn),a',':b','MarkerSize',20);
xlim([0.5 2.5]); xticks([1:2]); xticklabels({'acsf','nAChR'});
ylabel('latency to maximum (ms)'); ylim([0 200]); yticks([0:100:500]);
p = []; [~,p] = ttest2(a(:,1),a(:,2));
title(sprintf('DA to ACh peak - acsf:%d dhbe:%d ms \n latency ranksum: %1.4f | ttest: %1.4f',round(nanmean(a(:,1))),round(nanmean(a(:,2))),ranksum(a(:,1),a(:,2)),p)); 
axis('square');

subplot(2,2,2); hold on
a = v_diff([1 2 4],:); 
errorbar(nanmean(a)',SEM(a,1)','.k','MarkerSize',20);
plot([1.15; 1.85].*ones(2,nAn),a','.k','MarkerSize',20);
plot([1.15; 1.85].*ones(2,nAn),a',':b','MarkerSize',20);
xlim([0.5 2.5]); xticks([1:2]); xticklabels({'acsf','nAChR'});
ylabel('rDA1m peak amplitude (%dF/F)'); ylim([0 6]); yticks([-2:2:6]);
axis('square');
p = []; [~,p] = ttest2(a(:,1),a(:,2));
title(sprintf('DA to ACh peak - acsf:%1.2f dhbe:%1.2f dF/F \n max: ranksum: %1.4f | ttest: %1.4f',(nanmean(a(:,1))),(nanmean(a(:,2))),ranksum(a(:,1),a(:,2)),p)); 

movegui(gcf,'center');

% PLOT: AVERAGES for all groups
clr = {'k','r'};
% fig = figure; fig.Position([3 4]) = [1000 420];
subplot(2,2,3);
plot([0 0],[-6 6],'--k'); hold on; 
shadederrbar(time_seg, nanmean(da2ach_full{1},2), SEM(da2ach_full{1},2), clr{1}); 
shadederrbar(time_seg, nanmean(da2ach_full{2},2), SEM(da2ach_full{2},2), clr{2});
xlabel('Latency to ACh peak (s)'); xlim([-0.5 0.5]);
ylabel('rDA1m (%dF/F)'); ylim([-6 4]); yticks([-6:2:4])
axis('square')
title(sprintf('DA to ACh peak (n = %d mice)',nAn));

subplot(2,2,4);
plot([0 0],[-4 14],'--k'); hold on; 
shadederrbar(time_seg, nanmean(ach2ach_full{1},2), SEM(ach2ach_full{1},2), clr{1}); 
shadederrbar(time_seg, nanmean(ach2ach_full{2},2), SEM(ach2ach_full{2},2), clr{2});
xlabel('Latency to ACh peak (s)'); xlim([-0.5 0.5]);
ylabel('ACh3.0 (%dF/F)'); ylim([-4 10])
axis('square')
title(sprintf('ACh to ACh peak (n = %d mice)',nAn));
movegui(gcf,'center');

%% RELATIVE TO aCSF
nAn = 3;
fig = figure; fig.Position([3 4]) = [1000 900];
subplot(2,2,1); hold on
a = t_max([1 2 4],:)./nanmean(t_max([1 2 4],1));
plot([0.5 2.5],[1 1],'--k');
errorbar(abs(nanmean(a)'),abs(SEM(a,1)'),'.k','MarkerSize',20);
plot([1.15; 1.85].*ones(2,nAn),a','.k','MarkerSize',20);
plot([1.15; 1.85].*ones(2,nAn),a',':b','MarkerSize',20);
xlim([0.5 2.5]); xticks([1:2]); xticklabels({'aCSF','DHbE'});
ylabel('latency (max)'); ylim([0 2]); yticks([0:0.5:3]);
p = []; [~,p] = ttest(a(:,1),a(:,2));
title(sprintf('aCSF: %d ms || nAChR: %d ms \n latency signrank: %1.4f | ttest: %1.4f',round(nanmean(t_max([1 2 4],1))*1000),round(nanmean(t_max([1 2 4],2))*1000),signrank(a(:,1),a(:,2)),p)); 
axis('square');

subplot(2,2,2); hold on
a = v_diff([1 2 4],:)./nanmean(v_diff([1 2 4],1));
plot([0.5 2.5],[1 1],'--k');
errorbar(abs(nanmean(a)'),abs(SEM(a,1)'),'.k','MarkerSize',20);
plot([1.15; 1.85].*ones(2,nAn),a','.k','MarkerSize',20);
plot([1.15; 1.85].*ones(2,nAn),a',':b','MarkerSize',20);
xlim([0.5 2.5]); xticks([1:2]); xticklabels({'aCSF','DHbE'});
ylabel('difference in amplitude (max-min)'); ylim([0 2]); yticks([0:0.5:3]);
p = []; [~,p] = ttest(a(:,1),a(:,2));
title(sprintf('aCSF %1.2f dF/F || nAChR %1.2f dF/F \n maximum signrank: %1.4f | ttest: %1.4f',nanmean(v_diff([1 2 4],1)),nanmean(v_diff([1 2 4],2)),signrank(a(:,1),a(:,2)),p)); 
movegui(gcf,'center');

% PLOT: AVERAGES for all groups
clr = {'k','r'};
% fig = figure; fig.Position([3 4]) = [1000 420];
subplot(2,2,3);
plot([0 0],[-6 6],'--k'); hold on; 
shadederrbar(time_seg, nanmean(da2ach_full{1},2), SEM(da2ach_full{1},2), clr{1}); 
shadederrbar(time_seg, nanmean(da2ach_full{2},2), SEM(da2ach_full{2},2), clr{2});
xlabel('Latency to ACh peak (s)'); xlim([-0.5 0.5]); xticks([-0.5:0.25:0.5]);
ylabel('rDA1m (%dF/F)'); ylim([-4 4]); yticks([-6:2:4])
legend({'','aCSF','','DHbE'});
axis('square')
title(sprintf('DA to ACh peak (n = %d mice)',nAn));

subplot(2,2,4);
plot([0 0],[-4 14],'--k'); hold on; 
shadederrbar(time_seg, nanmean(ach2ach_full{1},2), SEM(ach2ach_full{1},2), clr{1}); 
shadederrbar(time_seg, nanmean(ach2ach_full{2},2), SEM(ach2ach_full{2},2), clr{2});
xlabel('Latency to ACh peak (s)'); xlim([-0.5 0.5]); xticks([-0.5:0.25:0.5]);
ylabel('ACh3.0 (%dF/F)'); ylim([-4 10])
axis('square')
title(sprintf('ACh to ACh peak (n = %d mice)',nAn));
movegui(gcf,'center');