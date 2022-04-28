%% (1) Load data
fPath = 'C:\Users\Anya\Desktop\IV_LOCAL\data_comb';
load(fullfile(fPath,'fullbeh_ACh_v2')); 
beh = behACh([14 16]); % Extract ACh recordings
load(fullfile(fPath,'fullcin_Feb21_ACh+DA+PF'));
sub = cinwt([20:25 28]);  % Extract CINs from recordings with ACh

%% CHECK
sm = 10;
nUnits = length([mat.n]); plm = floor(sqrt(nUnits)); pln = ceil(nUnits/plm);
% figure; for x = 1:nUnits; sp(x) = subplot(plm,pln,x); plot(time, movmean(align_rest(:,x),sm)); end
figure; for x = 1:nUnits; sp(x) = subplot(plm,pln,x); hold on
    shadederrbar(time, movmean(shuff50_rest(:,x),sm), movmean(shuff5_rest(:,x),sm), 'k');
    plot(time, movmean(alignDelta_rest(:,x),sm), 'b'); 
end
% figure; plot(time, movmean(nanmean(alignDelta_rest,2),sm));

%% (2) PLOT
fig = figure; fig.Position(3) = 1375;

subplot(1,3,1); hold on
sm = 10;
[m,i] = max(alignDelta_rest([1:50],:),[],1);
plot([0 0],[-0.2 0.2],'--k');
shadederrbar(time, movmean(nanmean(shuff50_rest,2),sm), movmean(nanmean(shuff5_rest,2),sm), 'k');
shadederrbar(time, movmean(nanmean(alignDelta_rest,2),sm), movmean(SEM(alignDelta_rest,2),sm), 'r'); 
ylabel('CIN Firing Rate (deltaFR)');
xlabel('Latency to ACh Peak (s)'); %grid on;
title(sprintf('max %1.1f +/- %1.1f || %d +/- %d ms',100*nanmean(m),100*SEM(m,2),round(1000*nanmean(time(i))),round(1000*SEM(time(i),1))));
axis('square');

above95 = []; below5 = [];
for y = 1:size(alignDelta_rest,2)
    a = alignDelta_rest(:,y); 
    a = a - nanmean(a([1:find(time == -0.51)],:));
    b = shuff95_rest(:,y);
    c = shuff50_rest(:,y) - (shuff95_rest(:,y) - shuff50_rest(:,y));  % mat(x).shuff5_rest(:,y);
    above95 = [above95, a > b]; %binary vector where CCG passed 95% confidence interval
    below5 = [below5, a < c]; %binary vector where CCG below 5% confidence interval
end
subplot(1,3,2); hold on
ds = 2;
a = 100*sum(above95,2)/size(above95,2); 
[prop_m,prop_t] = max(a); prop_t = time(prop_t); 
a = a(1:ds:end);
b = -100*sum(below5,2)/size(below5,2); b = b(1:ds:end);
bar(time(1:ds:end), a,'FaceColor','b','FaceAlpha',0.5);
bar(time(1:ds:end), b,'FaceColor','b','FaceAlpha',0.5);
xlabel('Lag (s)'); ylabel('Prop of Pair CCG > 95% CI (%)'); 
ylim([-100 100]); xlim([-1 1]);
title(sprintf('max %1.1f @ %d ms (n = %d)',prop_m, round(1000*prop_t), size(above95,2)))
axis('square');
movegui(gcf,'center');

%%
figure;
h = heatmap(alignDelta_rest');
h.Title = 'CIN to ACh peak';
h.XLabel = 'Latency to ACh Peak (s)'; h.YLabel = 'Unit';
h.GridVisible = 'off';
h.Colormap = jet; h.ColorLimits = [-0.5 1.5];