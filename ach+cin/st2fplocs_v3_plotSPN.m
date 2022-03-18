%% (1) Load data
fPath = 'C:\Users\Anya\Desktop\IV_LOCAL\data_comb';
load(fullfile(fPath,'fullbeh_ACh_v2')); 
beh = behACh([2 3 4 8 10 11 13]); % Extract ACh recordings
load(fullfile(fPath,'fullWT_2-10-21.mat'));
sub = WT(find([WT.label] == 1)); % Extract MSNs

%%
% run st2fplocs_v2 code

%%
fr_imm = [];
for x = 1:length(mat); fr_imm = [fr_imm; mat(x).fr(:,2)]; end

delta_new = alignDelta_rest;
tmp = delta_new(1,:) - shuff50_rest(1,:);
rmv = [find((abs(tmp) > 1.5)), find(fr_imm' > 4)];
delta_new(:,rmv) = nan;

delta_max = max(delta_new);
[delta_max,idx] = sort(delta_max);
delta_plot = delta_new(:,idx);

figure;
h = heatmap(delta_plot');
h.Title = 'MSN to ACh peak';
h.XLabel = 'Latency to ACh Peak (s)'; h.YLabel = 'Unit';
h.GridVisible = 'off';
h.Colormap = jet; h.ColorLimits = [-1 8];

%%
fig = figure; fig.Position(3) = 1000;
subplot(1,2,1); hold on
sm = 10;
% shuff50_adj = shuff50_rest(:,idx) - nanmean(shuff50_rest([1:50],idx)); % subtract baseline
shadederrbar(time, movmean(nanmean(shuff50_rest(:,idx),2),sm), movmean(nanmean(shuff5_rest(:,idx),2),sm), 'k');
shadederrbar(time, movmean(nanmean(delta_plot,2),sm), movmean(SEM(delta_plot,2),sm), 'b');
ylabel('pMSN firing rate change (%)'); ylim([-1 1])
title(sprintf('pMSN to ACh peak (n = %d units)', length(find(~isnan(delta_plot(1,:))))))
axis('square');

%
above95 = []; below5 = [];
for y = idx
    a = delta_new(:,y);
    b = shuff95_rest(:,y);
    c = shuff5_rest(:,y);
    above95 = [above95, a > b]; %binary vector where CCG passed 95% confidence interval
    below5 = [below5, a < c]; %binary vector where CCG below 5% confidence interval
end

% figure; hold on
subplot(1,2,2); hold on
ds = 2;
a = 100*sum(above95,2)/size(above95,2); 
[prop_m,prop_t] = max(a); prop_t = time(prop_t); 
a = a(1:ds:end);
b = -100*sum(below5,2)/size(below5,2); b = b(1:ds:end);
bar(time(1:ds:end), a,'FaceColor','b','FaceAlpha',0.5);
bar(time(1:ds:end), b,'FaceColor','b','FaceAlpha',0.5);
xlabel('Lag (s)'); ylabel('Prop of Pair CCG > 95% CI (%)'); 
ylim([-100 100]); yticks([-100:50:100]); xlim([-1 1]);
title(sprintf('max %1.1f @ %d ms',prop_m, round(1000*prop_t)))
axis('square');
movegui(gcf,'center');