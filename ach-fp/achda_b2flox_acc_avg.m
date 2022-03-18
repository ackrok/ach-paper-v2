a = [1:length(b2flox)]; rmv = [12 16 18 20 21 26]; a(rmv) = [];
a = [5 10 7 15 22 24 27];
beh = b2flox(a); % Extract recordings with reward
[align_full, time] = plot_fp2event(beh,[-4 1],0); % Align photometry to events

% a = [1:length(modAChDA)]; rmv = [3 5 17:25 27 33 30 36 41 42]; a(rmv) = [];
% beh = modAChDA(a); % Extract recordings with reward
% [align_full, time] = plot_fp2event(beh,[-4 1],0); % Align photometry to events

%%
align_adj = align_full;
for x = 1:length(beh); for y = 1:2
    align_adj{x,y} = align_adj{x,y} - nanmean(beh(x).FP{y});
    align_adj{x,y} = align_adj{x,y} - nanmean(align_adj{x,y}(find(time == -4):find(time == -0.5),:));
    end; end

% n = X mice
tmp = {}; for x = 1:length(beh); tmp{x} = strtok(beh(x).rec,'-'); end
uni = unique(tmp); nAn = length(uni);
align_an = cell(length(uni),1); 
align_an_avg = cell(2,1);
for x = 1:nAn; for y = 1:2
    ii = find(strcmp(tmp,uni{x}));
    align_an{x,y} = [align_adj{ii,y}];
    align_an_avg{y}(:,x) = nanmean(align_an{x,y},2);
    
    end; end
align_an_norm = cell(2,1);
for y = 1:2; align_an_norm{y} = normalize(align_an_avg{y},1,'range'); end

%% Population AVERAGE
fig = figure; fig.Position(3) = 1000;
subplot(1,2,1); hold on
plot([0 0],[-3 5],'--k');
shadederrbar(time, nanmean(align_an_avg{1},2), SEM(align_an_avg{1},2), 'g');
shadederrbar(time, nanmean(align_an_avg{2},2), SEM(align_an_avg{2},2), 'm');
ylabel('fluorescence'); % ylim([-2 8]); yticks([-2:2:8]);
xlabel('latency to acc (s)'); xlim([-1 1]); xticks([-1:0.5:1]);
title(sprintf('(n = %d mice)',nAn));
axis('square');

subplot(1,2,2); hold on
plot([0 0],[0 1],'--k');
shadederrbar(time, nanmean(align_an_norm{1},2), SEM(align_an_norm{1},2), 'g');
shadederrbar(time, nanmean(align_an_norm{2},2), SEM(align_an_norm{2},2), 'm');
ylabel('fluorescence'); yticks([0:0.5:1]);
xlabel('latency to acc (s)'); xlim([-0.5 0.5]); xticks([-1:0.5:1]);
title(sprintf('(n = %d mice)',nAn));
axis('square');
movegui(fig,'center');

%%
hold on
plot(time, nanmean(align_ach_avg,2), 'k');
plot(time, nanmean(align_da_avg,2), 'k');

%% NORM AVERAGE
