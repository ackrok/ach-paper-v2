%% ACh leads DA
beh = modAChDA; beh(42) = [];

x = 34;

fp_ach = beh(x).FP{1}; Fs = beh(x).Fs; % extract photometry signal from structure
fp_mu = nanmean(fp_ach); % mean of entire photometry signal
fp_ach = fp_ach - fp_mu; % SUBTRACT baseline% from fp%, now centered on BASELINE
idx_rew = extractEventST([1:length(fp_ach)]', floor(beh(x).reward)-100, floor(beh(x).reward)+100, 1); % identify sample during reward
idx_imm = extractEventST([1:length(fp_ach)]', beh(x).onRest, beh(x).offRest, 1); % identify sample during locomotion
idx_imm_nonRew = idx_imm(~ismember(idx_imm, idx_rew)); % exclude reward, include rest
fp_ach_imm = fp_ach(idx_imm_nonRew); 
fp_da = beh(x).FP{2}; % dopamine
fp_da = fp_da - nanmean(fp_da);
fp_da_imm = fp_da(idx_imm_nonRew);

figure; hold on; 
plot(beh(x).time, fp_ach, 'g'); 
plot(beh(x).time, fp_da, 'm');
stem(beh(x).reward./Fs, 10.*ones(length(beh(x).reward),1), 'b');
lick = beh(x).lick(:)./Fs;
a = [diff(lick.*1000) > 50]; % Identify licks that are <50ms after previous lick
lick_sub = lick; lick_sub(1) = [];
lick = [lick(1); lick_sub(a)];
plot([lick, lick]',[10;11].*ones(2,length(lick)),'-k');

ylabel('Photometry (%dF/F)'); % ylim([-20 20])

yyaxis right; 
acc = getAcc(beh(x).vel); 
plot(beh(x).time, acc, 'k');
ylabel('Acceleration (cm/s^2)'); ylim([-2 6]);
title(sprintf('%s',beh(x).rec))