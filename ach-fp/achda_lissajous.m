x = 15;
 
fp_1 = beh(x).FP{1}; Fs = beh(x).Fs; % extract photometry signal from structure
fp_mu = nanmean(fp_1); % mean of entire photometry signal
fp_1 = fp_1 - fp_mu; % SUBTRACT baseline% from fp%, now centered on BASELINE
fp_2 = beh(x).FP{2}; % dopamine
fp_2 = fp_2 - nanmean(fp_2);
f = [0.5 4];
fp_1_filt = filterFP(fp_1,Fs,f,10,'bandpass');
fp_2_filt = filterFP(fp_2,Fs,f,10,'bandpass');

win = [277 279]; % x = 30, immobility
win = [505 507]; % x = 30, locomotion
win = [1441.5 1442.6]; % x = 30, locomotion
win = [785 786]; % x = 15, reward
win = win.*Fs; nSamp = win(2) - win(1) + 1;

win_fp = [];
win_fp(:,1) = fp_1_filt(win(1):win(2));
win_fp(:,2) = fp_2_filt(win(1):win(2));

win_fp_norm = [];
win_fp_norm(:,1) = normalize(win_fp(:,1),'range');
win_fp_norm(:,2) = normalize(win_fp(:,2),'range');
win_fp_norm = win_fp_norm - nanmean(win_fp_norm);
%%
figure; hold on
plot(beh(x).time(win(1):win(2)), win_fp_norm);
plot(beh(x).time, [fp_1_filt, fp_2_filt]);
plot(beh(x).time, getAcc(beh(x).vel), 'k');
%%
fig = figure; fig.Position(3) = 1000;

clr = jet(nSamp);
subplot(2,2,1); hold on
plot(win_fp_norm(:,1),'k', 'Color', [0 0 0 0.2]);
for z = 1:nSamp
    plot(z,win_fp_norm(z,1),'.','MarkerSize',20,'Color',clr(z,:)); end
plot(win_fp_norm(:,2),'k'); 
ylabel('ACh intensity'); xlabel('time (s)'); xlim([0 nSamp]);
title(sprintf('%s: win=[%d %d] -- ACh',strtok(beh(x).rec,'R'),win(1)/Fs,win(2)/Fs));

subplot(2,2,3); hold on
plot(win_fp_norm(:,2),'k', 'Color', [0 0 0 0.2]);
for z = 1:nSamp
    plot(z,win_fp_norm(z,2),'.','MarkerSize',20,'Color',clr(z,:)); end
plot(win_fp_norm(:,1),'k'); 
ylabel('DA intensity'); xlabel('time (s)'); xlim([0 nSamp]);

subplot(2,2,[2 4]); hold on
plot(win_fp_norm(:,1), win_fp_norm(:,2), 'Color', [0 0 0 0.2]);
for z = 1:nSamp
    plot(win_fp_norm(z,1), win_fp_norm(z,2),'.','MarkerSize',20,'Color',clr(z,:)); end
xlabel('ACh intensity'); ylabel('DA intensity');