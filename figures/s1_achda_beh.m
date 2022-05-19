%% Extract recordings for BEHAVIOR
good_beh = [1:2 4:16 18:19 21 24:25 26 28:31 38:39 44 46:47];
beh = modAChDA(good_beh);
[align_full, time, ev] = plot_fp2event(beh,[-2 2],0); % Align photometry to events
Fs = 50; 

%% Behavior process
align_vel = cell(length(beh),1);
for x = 1:length(beh)
    ev_1 = ev{x}(~isnan(ev{x}));
    [sta_vel, t] = getSTA(beh(x).vel, ev_1, Fs, [-1 2]);
    align_vel{x} = sta_vel;
end

tmp = {}; for x = 1:length(beh); tmp{x} = strtok(beh(x).rec,'-'); end
uni = unique(tmp); nAn = length(uni);
vel_an = cell(length(uni),1); vel_an_avg = [];
for x = 1:nAn
    ii = find(strcmp(tmp,uni{x}));
    vel_an{x} = [align_vel{ii,1}];
    vel_an_avg(:,x) = nanmean([align_vel{ii,1}],2);
end

%% population AVERAGE
fig = figure; fig.Position(3) = 1000;
subplot(1,2,1); sm = 10;
shadederrbar(t, movmean(nanmean(vel_an_avg,2),sm), movmean(SEM(vel_an_avg,2),sm), 'k');
xlabel('Latency to onset (s)'); xlim([-1 2]); xticks([-1:0.5:2]);
ylabel('Velocity (cm/s)'); % ylim([0 1]); yticks([0:0.2:1]);
title(sprintf('Velocity Average (n = %d mice)',size(vel_an_avg,2)));

% LICK population HEATMAP
subplot(1,2,2);
% m = min(lick_an_avg);
% [a,b] = sort(m);
b = randperm(13,13);
h = heatmap(movmean(vel_an_avg(:,b)',sm,2));
h.Title = 'Velocity Average';
h.XLabel = 'Latency to onset (s)'; h.YLabel = 'Animal';
h.GridVisible = 'off';
h.Colormap = gray; % caxis([0 1.5])

%% SUBPLOTS
figure;
for x = 1:13
    sp(x) = subplot(4,4,x);
    a = movmean(vel_an{x},sm);
    plot(t, a, 'Color', [0 0 0 0.1]);
    shadederrbar(t, nanmean(a,2), SEM(a,2), 'b');
    % plot(t, nanmean(vel_an{x},2), 'k');
    title(sprintf('%s',uni{x}));
end
linkaxes(sp,'y');

%% COMBO PLOT
x = 3; % CHANGE
sm = 10;

fig = figure; fig.Position(3) = 1000;
subplot(1,2,1); hold on
a = movmean(vel_an{x},sm);
rmv = [find(a(t==1,:) < 0), find(abs(a(t==-0.24,:)) > 0.25)];
a(:,rmv) = [];
plot(t, a, 'Color', [0 0 0 0.1]);
shadederrbar(t, nanmean(a,2), SEM(a,2), 'b');
plot([0 0],[-1 20],'--k');
xlabel('Latency to onset (s)'); xlim([-1 2]); xticks([-1:0.5:2]);
ylabel('Velocity (cm/s)'); ylim([-1 20]); yticks([0:5:20]);
title(sprintf('%s: %d bouts',uni{x},size(a,2)));

subplot(1,2,2); hold on; 
a = movmean(vel_an_avg,sm);
% a(t < 0,:) = 0;
plot(t, a, 'Color', [0 0 0 0.1]);
shadederrbar(t, nanmean(a,2), SEM(a,2), 'b');
plot([0 0],[-1 12],'--k');
xlabel('Latency to onset (s)'); xlim([-1 2]); xticks([-1:0.5:2]);
ylabel('Velocity (cm/s)'); ylim([-1 12]); yticks([0:2:12]);
axis('square');