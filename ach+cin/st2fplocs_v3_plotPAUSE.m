%% (1) Load data
fPath = 'C:\Users\Anya\Desktop\IV_LOCAL\data_comb';
load(fullfile(fPath,'fullbeh_ACh_v2')); 
beh = behACh([2 3 4 8 10 11 13]); % Extract ACh recordings
load(fullfile(fPath,'fullcin_Feb21_ACh+DA+PF'));
sub = cinwt;  % Extract CINs from recordings with ACh

%% (2) PLOT
fig = figure; fig.Position(3) = 1375;

fr_imm = []; for x = 1:length(mat); fr_imm = [fr_imm; mat(x).fr(:,2)]; end

subplot(1,3,1); hold on
y = [13 14];
sm = 10;
plot([0 0],[3 7],'--k');
for z = 1:2
    % mid = (shuff50_rest(:,y(z)).*fr_imm(y(z))) + fr_imm(y(z));
    % high = (shuff95_rest(:,y(z)).*fr_imm(y(z))) + fr_imm(y(z));
    % low = (shuff5_rest(:,y(z)).*fr_imm(y(z))) + fr_imm(y(z));
    % shadederrbar(time, movmean(mid,sm), movmean(high-mid,sm), clr{z}); 
    plot(time, movmean(align_rest(:,y(z)),sm), 'k');
end
ylabel('Firing Rate (Hz)'); yticks([3:7]);
xlabel('Latency to ACh Pause (s)'); xlim([-1 1]); xticks([-1:0.5:1]);
title(sprintf('%s','IV066-rec03'));
axis('square');

subplot(1,3,2); hold on
a = alignDelta_rest - nanmean(alignDelta_rest([1:find(time == -0.51)],:)); % subtract baseline
b = shuff50_rest - nanmean(shuff50_rest([1:find(time == -0.51)],:)); % subtract baseline
shadederrbar(time, movmean(nanmean(b,2),sm), movmean(nanmean(shuff5_rest,2),sm), 'k');
shadederrbar(time, movmean(nanmean(a,2),5), movmean(SEM(a,2),5), 'r'); 
ylabel('CIN Firing Rate (deltaFR)'); ylim([-0.2 0.1]);
xlabel('Latency to ACh Pause (s)'); xlim([-1 1]); xticks([-1:0.5:1]);
title(sprintf('min = %1.2f',100*min(nanmean(a([1:end],:),2))));
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
subplot(1,3,3); hold on
ds = 2;
a = 100*sum(above95,2)/size(above95,2); a = a(1:ds:end);
b = -100*sum(below5,2)/size(below5,2); 
[prop_m,prop_t] = min(b); prop_t = time(prop_t); 
b = b(1:ds:end);
bar(time(1:ds:end), a,'FaceColor','b','FaceAlpha',0.5);
bar(time(1:ds:end), b,'FaceColor','b','FaceAlpha',0.5);
xlabel('Lag (s)'); ylabel('Prop of Pair CCG > 95% CI (%)'); 
ylim([-100 100]); xlim([-1 1]);
title(sprintf('max %1.1f @ %d ms (n = %d pairs)',prop_m, round(1000*prop_t), size(above95,2)))
axis('square');
movegui(gcf,'center');

% subplot(1,3,3); hold on
% sm = 1;
% [m,i] = max(alignDelta_rest([1:50],:),[],1);
% plot([0 0],[-0.2 0.6],'--k');
% shadederrbar(time, movmean(nanmean(shuff50_rest,2),sm), movmean(nanmean(shuff5_rest,2),sm), 'k');
% shadederrbar(time, movmean(nanmean(alignDelta_rest,2),sm), movmean(SEM(alignDelta_rest,2),sm), 'r'); 
% ylabel('CIN Firing Rate (deltaFR)');
% xlabel('Latency to ACh Peak (s)'); %grid on;
% title(sprintf('max %1.1f +/- %1.1f || %d +/- %d ms',100*nanmean(m),100*SEM(m,2),round(1000*nanmean(time(i))),round(1000*SEM(time(i),1))));
% axis('square');

%%
a = min(alignDelta_rest);
[~,b] = sort(a);
c = alignDelta_rest(:,flip(b));

figure;
h = heatmap(c');
h.Title = 'CIN to ACh pause';
h.XLabel = 'Latency to ACh Pause (s)'; h.YLabel = 'Unit';
h.GridVisible = 'off';
h.Colormap = jet; h.ColorLimits = [-0.4 0.4];