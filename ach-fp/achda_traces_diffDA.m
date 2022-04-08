x = 30;

fp_1 = beh(x).FP{1}; Fs = beh(x).Fs; % extract photometry signal from structure
fp_mu = nanmean(fp_1); % mean of entire photometry signal
fp_1 = fp_1 - fp_mu; % SUBTRACT baseline% from fp%, now centered on BASELINE
fp_2 = beh(x).FP{2}; % dopamine
fp_2 = fp_2 - nanmean(fp_2);
f = [0.5 4];
fp_1_filt = filterFP(fp_1,Fs,f,10,'bandpass');
fp_2_filt = filterFP(fp_2,Fs,f,10,'bandpass');

% fp_2_sub = [fp_2_filt(1); diff(fp_2_filt)]; % first derivative
tmp = hilbert(fp_2_filt);
fp_2_sub = complex(imag(tmp), -real(tmp)); % phase advance by 90 degrees by multiplying the signal by -1*i

figure; hold on
plot(beh(x).time, fp_1_filt, 'g', 'DisplayName', 'ACh');
plot(beh(x).time, fp_2_sub, 'm', 'DisplayName', 'diff(DA)');
ylabel('Photometry (%dF/F)'); ylim([-20 20])
yyaxis right; 
acc = getAcc(beh(x).vel); 
plot(beh(x).time, acc, 'k', 'DisplayName', 'acc');
ylabel('Acceleration (cm/s^2)'); ylim([-2 6]);
xlabel('Time (s)');
title(sprintf('%s',beh(x).rec)); legend
xlim([271 281]);

%%
x = 30;

fp_1 = beh(x).FP{1}; Fs = beh(x).Fs; % extract photometry signal from structure
fp_mu = nanmean(fp_1); % mean of entire photometry signal
fp_1 = fp_1 - fp_mu; % SUBTRACT baseline% from fp%, now centered on BASELINE
fp_2 = beh(x).FP{2}; % dopamine
fp_2 = fp_2 - nanmean(fp_2);
f = [0.5 4];
fp_1_filt = filterFP(fp_1,Fs,f,10,'bandpass');
fp_2_filt = filterFP(fp_2,Fs,f,10,'bandpass');

fp_2_sub = [fp_2_filt(1); diff(fp_2_filt)]; % first derivative

figure; hold on
plot(beh(x).time, fp_1_filt, 'g', 'DisplayName', 'ACh');
plot(beh(x).time, fp_2_filt, 'm', 'DisplayName', 'DA');
ylabel('Photometry (%dF/F)'); ylim([-20 20])
yyaxis right; 
acc = getAcc(beh(x).vel); 
plot(beh(x).time, acc, 'k', 'DisplayName', 'acc');
ylabel('Acceleration (cm/s^2)'); ylim([-2 6]);
xlabel('Time (s)');
title(sprintf('%s',beh(x).rec)); legend
xlim([271 281]);