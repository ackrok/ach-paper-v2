%% Input Variables
load('C:\Users\Anya\Desktop\IV_LOCAL\data_comb\fullcin_2022-Feb_WT');
beh = behWT([11:22, 32:45]); % DLS AP orientation recordings only
sub = cinWT;

%% Extract Peaks
for x = 1:length(beh)
    acc = getAcc(beh(x).vel); % Extract acceleration signal
    [~,locs] = findpeaks(acc,'MinPeakProminence',0.5,'MinPeakDistance',0.5); % Location of peaks, using findpeaks function
    time = [1/50:1/50:length(acc)/50]';
    beh(x).acc_locs = time(locs); % Convert peak locations to seconds
end; clc

%% Align CIN spikes to acceleration peaks
Fs = 50;
nShuff = 10; 
bin = 0.02; window = [-1 1]; %CHANGE: window for PETH 

mat = struct; % Initialize output structure
h = waitbar(0, 'PETH: CIN spikes to acceleration peaks');
for x = 1:length(beh)
    idx = find(strcmp({sub.rec},beh(x).rec)); % Find matching units from beh recording
    if isempty(idx); continue; end
    time = [1/50:1/50:length(beh(x).vel)/50]';
    st = {sub(idx).st}; % Unit spike times
    events = beh(x).acc_locs; % Acceleration peaks
    peth = getClusterPETH(st, events, bin, window); % PETH: CIN spikes to acc peaks
    
    mat(x).rec = beh(x).rec;
    mat(x).n = idx;
    mat(x).align = peth.fr;
    
    fr_mov = []; delta = [];
    tmp5 = []; tmp50 = []; tmp95 = []; 
    for y = 1:length(st)
        fr_mov(y) = 1/mean(diff(extractEventST(st{y},beh(x).on,beh(x).off,0))); % Firing rate during locomotion
        delta(:,y) = (peth.fr(:,y) - fr_mov(y))./fr_mov(y);
        % st_new = shuffleST(st{y},nShuff);
        % st_new = shiftST(st{y},nShuff,1000/nShuff); 
        st_new = poissonSpikeGen(sub(idx(y)).fr, time(end), nShuff);
        peth_shuff = getClusterPETH(st_new, events, bin, window);
        prc = prctile(peth_shuff.fr,[5 50 95],2); %5th, 50th, 95th percentile of shuffled PETH
        prc = (prc - fr_mov(y))./ fr_mov(y); % Convert to delta
        tmp5(:,y) = prc(:,1); tmp50(:,y) = prc(:,2); tmp95(:,y) = prc(:,3);
    end
    mat(x).fr = fr_mov;
    mat(x).delta = delta;
    mat(x).prc5 = tmp5;
    mat(x).prc50 = tmp50;
    mat(x).prc95 = tmp95;
    waitbar(x/length(beh),h);
end; fprintf('Done! \n'); close(h);
time = peth.time;

%% EXTRACT
nUnits = length([mat.n]);
delta_full = reshape([mat.delta],[length(peth.time),nUnits]); 
shuff5 = reshape([mat.prc5],[length(peth.time),nUnits]); 
shuff50 = reshape([mat.prc50],[length(peth.time),nUnits]); 
shuff95 = reshape([mat.prc95],[length(peth.time),nUnits]); 

%% PLOT AVERAGE
% figure; hold on
fig = figure; fig.Position(3) = 1375;
subplot(1,3,1); hold on
sm = 2;
plot([0 0],[-0.5 0.5],'--k');
% shadederrbar(time, nanmean(movmean(shuff50,sm,1),2), nanmean(movmean(shuff95-shuff50,sm,1),2), 'k'); 
shadederrbar(time, nanmean(movmean(delta_full,sm,1),2), SEM(movmean(delta_full,sm,1),2), 'b'); 
xlabel('Latency (s)'); xlim([-1 1]); xticks([-1:0.5:1]);
ylabel('delta firing rate'); ylim([-0.5 0.5]); yticks([-0.5:0.25:0.5]);
title(sprintf('st to acc (n = %d units)',nUnits)); axis('square');

%% SUBPLOT delta
figure;
for x = 1:length(mat)
    sp(x) = subplot(6,6,x);
    plot(time, mat(x).delta)
    title(sprintf('%s',mat(x).rec));
end
linkaxes(sp,'y');

%% HEATMAP
a = max(delta_full); [~,b] = sort(a); c = delta_full(:,b);
% figure; 
h = heatmap(c');
h.Title = 'CIN to acc peak';
h.XLabel = 'Latency to acc (s)'; h.YLabel = 'Unit';
h.GridVisible = 'off';
h.Colormap = jet; h.ColorLimits = [-1 2];

%% PLOT proportion > 95% CI
above95 = []; below5 = [];
for x = 1:nUnits
    a = delta_full(:,x); b = shuff95(:,x); c = shuff5(:,x);
    above95 = [above95, a > b]; %binary vector where CCG passed 95% confidence interval
    below5 = [below5, a < c]; %binary vector where CCG below 5% confidence interval
end
subplot(1,3,2); hold on % figure; hold on
ds = 2;
a = 100*sum(above95,2)/size(above95,2); 
[prop_m,prop_t] = max(a); prop_t = time(prop_t); 
a = a(1:ds:end);
b = -100*sum(below5,2)/size(below5,2); b = b(1:ds:end);
bar(time(1:ds:end), a,'FaceColor','b','FaceAlpha',0.5);
bar(time(1:ds:end), b,'FaceColor','b','FaceAlpha',0.5);
xlabel('Lag (s)'); ylabel('Prop of Pair CCG > 95% CI (%)'); 
ylim([-100 100]); xlim([-1 1]);
title(sprintf('max %1.1f @ %d ms (n = %d units)',prop_m, round(1000*prop_t), size(above95,2)))
axis('square');

