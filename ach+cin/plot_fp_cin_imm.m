%% (1) Load data
fPath = 'C:\Users\Anya\Desktop\IV_LOCAL\data_comb';
load(fullfile(fPath,'fullbeh_ACh_v2')); 
beh = behACh([2 3 4 8 10 11 13]); % Extract ACh recordings
load(fullfile(fPath,'fullcin_Feb21_ACh+DA+PF'));
sub = cinwt;  % Extract CINs from recordings with ACh

%% ACh + CIN
x = 4; % recording
z = 2; % immobility

idx = strcmp({sub.rec},beh(x).rec);
st = {sub(idx).st}; % Extract spike times of units from this recording
fp = beh(x).FP{1};
fp = fp - nanmean(fp); % subtract mean of baseline recording to center around 0%
thres_peak = 2; % CHANGE -- threshold for peaks/pauses to exceed, after baseline
thres_pause = -1;
width_peak = 5; % CHANGE -- require that peaks/pauses are at least X samples long (1 sample = 20ms)
locs = findIdxPause(fp,'peak',thres_peak,width_peak);
locs_imm = extractEventST(locs, beh(x).onRest, beh(x).offRest, 1); % Event times during rest
locs_imm = locs_imm./beh(x).Fs; % convert to seconds
idx_imm = extractEventST([1:length(fp)]', beh(x).onRest, beh(x).offRest, 1); % Event times during rest
bin_imm = zeros(length(fp),1); bin_imm(idx_imm) = 1;

figure; hold on
area(beh(x).time, 10.*bin_imm,'FaceAlpha',0.2,'FaceColor','r','EdgeAlpha',0);
plot(beh(x).time, movmean(fp,5), 'Color', [0.05 0.75 0.45]);
plot(beh(x).time, -2+getAcc(beh(x).vel), 'k');
% stem(locs_imm, 5.*ones(length(locs_imm),1), 'k');
for y = 1:length(st)
    plot(st{y},y.*ones(length(st{y}),1),'.b','MarkerSize',20);
end
title(sprintf('%s',beh(x).rec),'Interpreter','none');
ylabel('ACh3.0 fluorescence'); % ylim([-4 8]); yticks([-4:8]); 
xticks([1:0.5:6000]);
% xlim([2165 2173.5]); % x = 2
% xlim([4908.5 4914.5]); % x = 11
% xlim([5025.5 5032.5]); % x = 11
% xlim([1989 1995]); x = 4
xlim([200 300]);

%% DA + CIN
x = 9; % recording
z = 2; % immobility

idx = find(strcmp({sub.rec},beh(x).rec));
st = {sub(idx).st}; % Extract spike times of units from this recording
fp = beh(x).FP{1};
fp = fp - nanmean(fp); % subtract mean of baseline recording to center around 0%
thres_peak = 2; % CHANGE -- threshold for peaks/pauses to exceed, after baseline
thres_pause = -1;
width_peak = 5; % CHANGE -- require that peaks/pauses are at least X samples long (1 sample = 20ms)
locs = findIdxPause(fp,'peak',thres_peak,width_peak);
locs_imm = extractEventST(locs, beh(x).onRest, beh(x).offRest, 1); % Event times during rest
locs_imm = locs_imm./beh(x).Fs; % convert to seconds
idx_imm = extractEventST([1:length(fp)]', beh(x).onRest, beh(x).offRest, 1); % Event times during rest
bin_imm = zeros(length(fp),1); bin_imm(idx_imm) = 1;

figure; hold on
area(beh(x).time, 10.*bin_imm,'FaceAlpha',0.2,'FaceColor','r','EdgeAlpha',0);
plot(beh(x).time, movmean(fp,5), 'Color', [0.05 0.75 0.45]);
plot(beh(x).time, -2+getAcc(beh(x).vel), 'k');
% stem(locs_imm, 5.*ones(length(locs_imm),1), 'k');
for y = 1:length(st)
    plot(st{y},y.*ones(length(st{y}),1),'.b','MarkerSize',20);
end
title(sprintf('%s',beh(x).rec),'Interpreter','none');
ylabel('DA4.4 fluorescence'); legend({'rest','DA4.4','acceleration','pCIN','pCIN','pCIN'});
ylim([-10 10]); yticks([-10:5:10]); xticks([1:0.5:2500]);
xlim([314 324]);