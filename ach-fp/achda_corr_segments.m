beh = modAChDA; beh(42) = [];

x = 41;
fp_mat = beh(x).FP{1} - nanmean(beh(x).FP{1}); % extract ACh photometry data
fp_mat(:,2) = beh(x).FP{2} - nanmean(beh(x).FP{2}); % extract DA photometry data

a = find(beh(x).time == [676]); % seg1: more 180deg out of phase
a(2) = find(beh(x).time == [681.5]);
seg1 = fp_mat([a(1):a(2)],:);
seg1 = seg1 - nanmean(seg1,1);

a = find(beh(x).time == [1581.5]); % seg2: appears slightly shifted
a(2) = find(beh(x).time == [1587.5]);
seg2 = fp_mat([a(1):a(2)],:);
seg2 = seg2 - nanmean(seg2,1);

Fs = 50;
corr_tmp = [];
[corr_tmp(:,1), lags] = xcorr(seg1(:,1), seg1(:,2), 2*Fs, 'coeff'); % cross-correlation
[corr_tmp(:,2), lags] = xcorr(seg2(:,1), seg2(:,2), 2*Fs, 'coeff'); % cross-correlation

fig = figure; fig.Position(3) = 1000;
subplot(2,2,1); hold on
plot(seg1(:,1), 'g'); plot(seg1(:,2), 'm');
xlabel('samples'); xlim([0 300]); ylabel('fluorescence (%dF/F)');
title('1');

subplot(2,2,3); hold on
plot(seg2(:,1), 'g'); plot(seg2(:,2), 'm');
xlabel('samples'); xlim([0 300]); ylabel('fluorescence (%dF/F)'); ylim([-10 15]);
title('2');

subplot(2,2,[2 4]);
plot(lags/Fs, corr_tmp);
xlabel('lags (s)'); ylabel('coeff'); ylim([-1 1]);
axis('square');
legend({'1','2'});