%% VISUALIZE
figure;
for x = 1:size(ccgDelta_rest,2)
    sp(x) = subplot(4,4,x); hold on
    sm = 5;
    shadederrbar(time, movmean(ccgDelta_50rest(:,x),sm), movmean(ccgDelta_95rest(:,x),sm), 'k'); hold on
    shadederrbar(time, movmean(ccgDelta_rest(:,x),sm), movmean(ccgDelta_rest(:,x),sm), 'b'); 
    % shadederrbar(time, movmean(ccgDelta_50mvmt(:,x),sm), movmean(ccgDelta_95mvmt(:,x),sm), 'g'); hold on
    % shadederrbar(time, movmean(ccgDelta_mvmt(:,x),sm), movmean(ccgDelta_mvmt(:,x),sm), 'b'); 
    xlabel('Lag (s)'); ylabel('CCG (deltaFR)'); xlim([-1 1]); grid on; 
    title(sprintf('%d)',x));
end
linkaxes(sp,'y');

%% HEATMAP
figure;
a = ccgDelta_rest;
b = a(time == 0,:); % find delta @ lag = 0
[c, ii] = sort(b); % sort in ascending order
h = heatmap(a([find(time == -0.5):find(time == 0.5)],ii)');
h.Title = 'REST CIN pair CCG';
h.XLabel = 'Lag from reference spike (s)'; 
h.YLabel = 'Pair Number';
h.Colormap = jet; h.GridVisible = 'off'; h.ColorLimits = [-0.2 1]; 

%%
figure;
a = ccgDelta_mvmt;
b = a(time == 0,:); % find delta @ lag = 0
% [c, ii] = sort(b); % sort in ascending order
h = heatmap(a([find(time == -0.5):find(time == 0.5)],ii)');
h.Title = 'MOV CIN pair CCG';
h.XLabel = 'Lag from reference spike (s)'; 
h.YLabel = 'Pair Number';
h.Colormap = jet; h.GridVisible = 'off'; h.ColorLimits = [-1 1.5];

%% (1) PLOT REST
fig = figure; fig.Position(3) = 1375;

y = [3 4 13]; clr = {'c','b'}; 
sm = 10;

subplot(1,3,1); hold on
plot([0 0],[3.5 5.5],'--k');
for z = 1:length(y)
%     mid = (shuff50_rest(:,y(z)).*fr_imm(y(z))) + fr_imm(y(z));
%     high = (shuff95_rest(:,y(z)).*fr_imm(y(z))) + fr_imm(y(z));
%     low = (shuff5_rest(:,y(z)).*fr_imm(y(z))) + fr_imm(y(z));
%     shadederrbar(time, movmean(mid,sm), movmean(high-mid,sm), clr{z}); 
    plot(time, movmean(ccg_rest(:,y(z)),sm), 'k');
end
ylabel('Firing Rate (Hz)'); 
xlabel('Latency to ref. spike (s)'); xlim([-0.5 0.5]); xticks([-1:0.25:1]);
% title(sprintf('%s','IV059-rec02'));
axis('square');

subplot(1,3,2); hold on
a = ccgDelta_rest - nanmean(ccgDelta_rest([1:50],:)); % subtract baseline
b = ccgDelta_50rest - nanmean(ccgDelta_50rest([1:50],:)); % subtract baseline
c = ccgDelta_95rest - nanmean(ccgDelta_95rest([1:50],:)); % subtract baseline
shadederrbar(time, movmean(nanmean(ccgDelta_50rest,2),sm), movmean(nanmean(ccgDelta_95rest,2),sm), 'k');
% shadederrbar(time, movmean(nanmean(b,2),sm), movmean(nanmean(c,2),sm), 'k');
shadederrbar(time, movmean(nanmean(a,2),1), movmean(SEM(a,2),1), 'b'); 
ylabel('CIN Firing Rate (deltaFR)'); yticks([-0.1:0.1:0.3]);
xlabel('Latency to ref. spike (s)'); xlim([-0.5 0.5]); xticks([-1:0.25:1]);
title(sprintf('REST max = %1.3f',max(nanmean(a([51:end],:),2))));
axis('square');

above95 = []; below5 = [];
for y = 1:size(ccgDelta_rest,2)
    a = ccgDelta_rest(:,y); 
    b = ccgDelta_95rest(:,y);
    c = ccgDelta_50rest(:,y) - (ccgDelta_95rest(:,y) - ccgDelta_50rest(:,y));  % mat(x).shuff5_rest(:,y);
    above95 = [above95, a > b]; %binary vector where CCG passed 95% confidence interval
    below5 = [below5, a < c]; %binary vector where CCG below 5% confidence interval
end
subplot(1,3,3); hold on
ds = 2;
a = 100*sum(above95,2)/size(above95,2); 
[prop_m,prop_t] = max(a); prop_t = time(prop_t); 
a = a(1:ds:end);
b = -100*sum(below5,2)/size(below5,2); b = b(1:ds:end);
bar(time(1:ds:end), a,'FaceColor','b','FaceAlpha',0.5);
bar(time(1:ds:end), b,'FaceColor','b','FaceAlpha',0.5);
xlabel('Latency to ref. spike (s)'); ylabel('Prop of Pair CCG > 95% CI (%)'); 
ylim([-100 100]); xlim([-0.5 0.5]);
title(sprintf('max %1.1f @ %d ms (n = %d units)',prop_m, round(1000*prop_t), size(above95,2)))
axis('square');
movegui(gcf,'center');

%% (2) PLOT MOV
fig = figure; fig.Position(3) = 1375;

y = [3 4 13]; clr = {'c','b'}; 
sm = 10;

subplot(1,3,1); hold on
plot([0 0],[4 6],'--k');
for z = 1:length(y)
%     mid = (shuff50_rest(:,y(z)).*fr_imm(y(z))) + fr_imm(y(z));
%     high = (shuff95_rest(:,y(z)).*fr_imm(y(z))) + fr_imm(y(z));
%     low = (shuff5_rest(:,y(z)).*fr_imm(y(z))) + fr_imm(y(z));
%     shadederrbar(time, movmean(mid,sm), movmean(high-mid,sm), clr{z}); 
    plot(time, movmean(ccg_mvmt(:,y(z)),sm), 'k');
end
ylabel('Firing Rate (Hz)'); 
xlabel('Latency to ref. spike (s)'); xlim([-1 1]); xticks([-1:0.5:1]);
axis('square');

subplot(1,3,2); hold on
a = ccgDelta_mvmt - nanmean(ccgDelta_mvmt([1:50],:)); % subtract baseline
shadederrbar(time, movmean(nanmean(ccgDelta_50mvmt,2),sm), movmean(nanmean(ccgDelta_95mvmt,2),sm), 'k');
shadederrbar(time, movmean(nanmean(a,2),1), movmean(SEM(a,2),1), 'r'); 
ylabel('CIN Firing Rate (deltaFR)'); yticks([-0.5:0.1:0.5]);
xlabel('Latency to ref. spike (s)'); xlim([-1 1]); xticks([-1:0.5:1]);
title(sprintf('MOV max = %1.3f',max(nanmean(a([51:end],:),2))));
axis('square');

above95 = []; below5 = [];
for y = 1:size(ccgDelta_mvmt,2)
    a = ccgDelta_mvmt(:,y); 
    b = ccgDelta_95mvmt(:,y);
    c = ccgDelta_50mvmt(:,y) - (ccgDelta_95mvmt(:,y) - ccgDelta_50mvmt(:,y));  % mat(x).shuff5_rest(:,y);
    above95 = [above95, a > b]; %binary vector where CCG passed 95% confidence interval
    below5 = [below5, a < c]; %binary vector where CCG below 5% confidence interval
end
subplot(1,3,3); hold on
ds = 2;
a = 100*sum(above95,2)/size(above95,2); 
[prop_m,prop_t] = max(a); prop_t = time(prop_t); 
a = a(1:ds:end);
b = -100*sum(below5,2)/size(below5,2); b = b(1:ds:end);
bar(time(1:ds:end), a,'FaceColor','b','FaceAlpha',0.5);
bar(time(1:ds:end), b,'FaceColor','b','FaceAlpha',0.5);
xlabel('Latency to ref. spike (s)'); ylabel('Prop of Pair CCG > 95% CI (%)'); 
ylim([-100 100]); xlim([-1 1]);
title(sprintf('max %1.1f @ %d ms (n = %d units)',prop_m, round(1000*prop_t), size(above95,2)))
axis('square');
movegui(gcf,'center');