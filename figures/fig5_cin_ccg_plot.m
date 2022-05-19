%% LOAD DATA
cd('R:\homes\ack466\ACh paper\3_CIN_rest_synchrony');
% data_cin_ccg_rest+mov_v4: DLS pCINs
%   2022-05-18: updated 95% CI to be 2.5% to 97.5%
%
% data_cin_ccg_rest+mov_DMS_v4: DMS pCINs
%   2022-05-18: updated 95% CI to be 2.5% to 97.5%
%
% data_cin_ccg_rest+mov_6ohda
%   NOT YET UPDATED 95% CI
%
% data_cin_ccg_rest+mov_chat
%   NOT YET UPDATED 95% CI
%
% data_cin_ccg_rest+mov_gat
%   NOT YET UPDATED 95% CI
%

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

%% EXTRACT FROM MAT
time = [-2:0.01:2];
ccg_full = []; ccg_mvmt = []; ccg_rest = [];
ccgDelta = []; ccg95 = []; ccg50 = [];
ccgDelta_mvmt = []; ccgDelta_95mvmt = []; ccgDelta_50mvmt = [];
ccgDelta_rest = []; ccgDelta_95rest = []; ccgDelta_50rest = [];

for x = 1:length(mat)
    if isempty(mat(x).ccg_rest); continue; end
    ccg_full = [ccg_full, mat(x).ccg];
    ccg_mvmt = [ccg_mvmt, mat(x).ccg_mvmt]; 
    ccg_rest = [ccg_rest, mat(x).ccg_rest];
    for y = 1:length(mat(x).fr_rest)
        tmp = (mat(x).ccg(:,y) - mat(x).fr(y))./mat(x).fr(y); %deltaFR = (rate - mu(rate))/mu(rate)
        ccgDelta = [ccgDelta, tmp];
        tmp = (mat(x).ccg_rest(:,y) - mat(x).fr_rest(y))./mat(x).fr_rest(y); %deltaFR = (rate - mu(rate))/mu(rate)
        ccgDelta_rest = [ccgDelta_rest, tmp];
        tmp = (mat(x).ccg_mvmt(:,y) - mat(x).fr_mvmt(y))./mat(x).fr_mvmt(y); %deltaFR = (rate - mu(rate))/mu(rate)
        ccgDelta_mvmt = [ccgDelta_mvmt, tmp];
        
        tmp_95 = (mat(x).shuffPrc{y}(:,3) - mat(x).fr(y))./mat(x).fr(y); 
        tmp_50 = (mat(x).shuffPrc{y}(:,2) - mat(x).fr(y))./mat(x).fr(y);
        ccg95 = [ccg95, tmp_95];
        ccg50 = [ccg50, tmp_50];
        
        tmp_95 = (mat(x).shuffPrc_rest{y}(:,3) - mat(x).fr_rest(y))./mat(x).fr_rest(y); 
        tmp_50 = (mat(x).shuffPrc_rest{y}(:,2) - mat(x).fr_rest(y))./mat(x).fr_rest(y);
        ccgDelta_95rest = [ccgDelta_95rest, tmp_95];
        ccgDelta_50rest = [ccgDelta_50rest, tmp_50];
        
        tmp_95 = (mat(x).shuffPrc_mvmt{y}(:,3) - mat(x).fr_mvmt(y))./mat(x).fr_mvmt(y); 
        tmp_50 = (mat(x).shuffPrc_mvmt{y}(:,2) - mat(x).fr_mvmt(y))./mat(x).fr_mvmt(y);
        ccgDelta_95mvmt = [ccgDelta_95mvmt, tmp_95];
        ccgDelta_50mvmt = [ccgDelta_50mvmt, tmp_50];
    end
end

%% IMMOBILITY 
for mm = 1
fig = figure; fig.Position(3) = 1375;

% subplot(1,3,1);
% a = ccgDelta_rest(:,[1:end]);
% b = a(time == 0,:); % find delta @ lag = 0
% [c, ii] = sort(b); % sort in ascending order
% h2 = heatmap(a([find(time == -0.5):find(time == 0.5)],ii)');
% h2.Title = 'Cre(+) CCG imm';
% h2.XLabel = 'Lag from ref. spike (s)'; 
% h2.YLabel = 'Pair Number';
% h2.Colormap = jet; h2.GridVisible = 'off'; % h2.ColorLimits = [-0.4 1.2]; 

subplot(1,3,1);
a = ccgDelta_rest(:,[1:end]);
b = a(time == 0,:); % find delta @ lag = 0
[c, ii] = sort(b); % sort in ascending order
d = a([find(time == -1):find(time == 1)],ii)';
h = imagesc(time, [1:size(a,2)], d, [0 1]);
colorbar; colormap(jet(256));
title('CCG imm'); axis('square')
xlabel('Lag from ref. spike (s)'); ylabel('Pair Number');
h3 = h.Parent; h3.CLim = [-0.2 0.8];

subplot(1,3,2);
sm = 5;
shadederrbar(time, movmean(nanmean(ccgDelta_50rest(:,[1:end]),2),sm), movmean(nanmean(ccgDelta_95rest(:,[1:end]),2),sm), 'r'); hold on
shadederrbar(time, movmean(nanmean(ccgDelta_rest(:,[1:end]),2),sm), movmean(SEM(ccgDelta_rest(:,[1:end]),2),sm), 'r'); 
xlabel('Lag from ref. spike (s)'); xlim([-1 1]); xticks([-1:0.5:1]);
ylabel('pCIN firing rate change (%)'); % ylim([-0.1 0.3]); yticks([-0.1:0.1:0.3]); 
title(sprintf('max %1.3f (n = %d pairs)',max(movmean(nanmean(ccgDelta_rest(:,[1:end]),2),sm)),size(a,2)));
axis('square');

above95 = []; below5 = [];
for x = 1:length(mat)
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
b = -100*sum(below5,2)/size(below5,2); 
[prop_m2] = min(b);
b = b(1:ds:end);
bar(time(1:ds:end), a,'FaceColor','b','FaceAlpha',0.5);
bar(time(1:ds:end), b,'FaceColor','b','FaceAlpha',0.5);
xlabel('Lag from ref. spike (s)'); xlim([-1 1]); xticks([-1:0.5:1]);
ylabel('Prop of pairs > 95% CI (%)'); ylim([-50 100]); xlim([-1 1]);
title(sprintf('max %1.1f @ %d ms | min %1.1f',prop_m, round(1000*prop_t), prop_m2))
axis('square');
movegui(gcf,'center');
end

%% LOCOMOTION
for mm = 1
fig = figure; fig.Position(3) = 1375;

% subplot(1,3,1);
% a = ccgDelta_rest(:,[1:end]);
% b = a(time == 0,:); % find delta @ lag = 0
% [c, ii] = sort(b); % sort in ascending order
% h2 = heatmap(a([find(time == -0.5):find(time == 0.5)],ii)');
% h2.Title = 'Cre(+) CCG imm';
% h2.XLabel = 'Lag from ref. spike (s)'; 
% h2.YLabel = 'Pair Number';
% h2.Colormap = jet; h2.GridVisible = 'off'; % h2.ColorLimits = [-0.4 1.2]; 

subplot(1,3,1);
a = ccgDelta_mvmt(:,[1:end]);
b = a(time == 0,:); % find delta @ lag = 0
[c, ii] = sort(b); % sort in ascending order
d = a([find(time == -1):find(time == 1)],ii)';
h = imagesc(time, [1:size(a,2)], d, [0 1]);
colorbar; colormap(jet(256));
title('CCG locomotion'); axis('square')
xlabel('Lag from ref. spike (s)'); ylabel('Pair Number');
h2 = h.Parent; h2.CLim = [-0.25 1.5];

subplot(1,3,2);
sm = 5;
shadederrbar(time, movmean(nanmean(ccgDelta_50mvmt(:,[1:end]),2),sm), movmean(nanmean(ccgDelta_95mvmt(:,[1:end]),2),sm), 'r'); hold on
shadederrbar(time, movmean(nanmean(ccgDelta_mvmt(:,[1:end]),2),sm), movmean(SEM(ccgDelta_mvmt(:,[1:end]),2),sm), 'r'); 
xlabel('Lag from ref. spike (s)'); xlim([-1 1]); xticks([-1:0.5:1]);
ylabel('pCIN firing rate change (%)'); % ylim([-0.1 0.3]); yticks([-0.1:0.1:0.3]); 
title(sprintf('CCG locomotion (n = %d pairs)',size(a,2)));
axis('square');

above95 = []; below5 = [];
for x = 1:length(mat)
    for y = 1:size(mat(x).ccg_mvmt,2)
        a = mat(x).ccg_mvmt(:,y); b = mat(x).shuffPrc_mvmt{y};
        above95 = [above95, a > b(:,3)]; %binary vector where CCG passed 95% confidence interval
        below5 = [below5, a < b(:,1)]; %binary vector where CCG below 5% confidence interval
    end
end
subplot(1,3,3); hold on
ds = 2;
a = 100*sum(above95,2)/size(above95,2); 
[prop_m,prop_t] = max(a); prop_t = time(prop_t); 
a = a(1:ds:end);
b = -100*sum(below5,2)/size(below5,2); prop_m2 = min(b); b = b(1:ds:end);
bar(time(1:ds:end), a,'FaceColor','b','FaceAlpha',0.5);
bar(time(1:ds:end), b,'FaceColor','b','FaceAlpha',0.5);
xlabel('Lag from ref. spike (s)'); xlim([-1 1]); xticks([-1:0.5:1]);
ylabel('Prop of pairs > 95% CI (%)'); ylim([-50 100]); xlim([-1 1]);
title(sprintf('max %1.1f @ %d ms | min %1.1f',prop_m, round(1000*prop_t), prop_m2))
axis('square');
movegui(gcf,'center');
end
