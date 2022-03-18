%% FIRING RATE
fr_imm = [];
for x = 1:length(sub)
    ib = find(strcmp({beh.rec},sub(x).rec));
    diffFs = 50;
    fr_imm(x) = 1/mean(diff(extractEventST(sub(x).st,beh(ib).onRest/diffFs,beh(ib).offRest/diffFs,0)));
end

figure; hold on
a = {fr_imm([1:32]), fr_imm([33:51])};
violinplot(a);
errorbar(cellfun(@nanmean,a), cellfun(@nanstd,a)./sqrt(cellfun(@length,a)), '.k', 'MarkerSize', 20);
ylabel('pCIN firing rate'); ylim([0 10]); yticks([0:2:10]);
xticklabels({'cre (-)','cre (+)'});
[~,p] = ttest2(a{1},a{2});
title(sprintf('ranksum: %1.4f, ttest2: %1.4f',ranksum(a{1},a{2}),p));

%% PLOT 
fig = figure; fig.Position(3) = 1375;

subplot(1,3,1);
a = ccgDelta_rest(:,[1:40]);
b = a(time == 0,:); % find delta @ lag = 0
[c, ii] = sort(b); % sort in ascending order
h2 = heatmap(a([find(time == -0.5):find(time == 0.5)],ii)');
h2.Title = 'Cre(+) CCG imm';
h2.XLabel = 'Lag from ref. spike (s)'; 
h2.YLabel = 'Pair Number';
h2.Colormap = jet; h2.GridVisible = 'off'; % h2.ColorLimits = [-0.4 1.2]; 

subplot(1,3,2);
sm = 5;
shadederrbar(time, movmean(nanmean(ccgDelta_50rest(:,[1:40]),2),sm), movmean(nanmean(ccgDelta_95rest(:,[1:40]),2),sm), 'r'); hold on
shadederrbar(time, movmean(nanmean(ccgDelta_rest(:,[1:40]),2),sm), movmean(SEM(ccgDelta_rest(:,[1:40]),2),sm), 'r'); 
xlabel('Lag from ref. spike (s)'); xlim([-1 1]); xticks([-1:0.5:1]);
ylabel('pCIN firing rate change (%)'); % ylim([-0.1 0.3]); yticks([-0.1:0.1:0.3]); 
title(sprintf('CCG imm (n = %d pairs)',size(a,2)));
axis('square');

above95 = []; below5 = [];
for x = 1:3
    for y = 1:size(mat(x).ccg_rest,2)
        a = mat(x).ccg_rest(:,y); b = mat(x).shuffPrc_rest{y};
        above95 = [above95, a > b(:,3)]; %binary vector where CCG passed 95% confidence interval
        below5 = [below5, a < b(:,1)]; %binary vector where CCG below 5% confidence interval
    end
end
subplot(1,3,3); hold on
ds = 2;
a = 100*sum(above95,2)/size(above95,2); 
[prop_m,prop_t] = max(a); prop_t = time(prop_t); 
a = a(1:ds:end);
b = -100*sum(below5,2)/size(below5,2); b = b(1:ds:end);
bar(time(1:ds:end), a,'FaceColor','b','FaceAlpha',0.5);
bar(time(1:ds:end), b,'FaceColor','b','FaceAlpha',0.5);
xlabel('Lag from ref. spike (s)'); xlim([-1 1]); xticks([-1:0.5:1]);
ylabel('Prop of pairs > 95% CI (%)'); ylim([-100 100]); xlim([-1 1]);
title(sprintf('max %1.1f @ %d ms (n = %d units)',prop_m, round(1000*prop_t), size(above95,2)))
axis('square');
movegui(gcf,'center');