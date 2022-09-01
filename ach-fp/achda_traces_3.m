%% ACh DA - figure 1 AK189
% 5 sec segments from (1) reward (2) locomotion (3) immobility
% 220506: (1) x = 43; xlim([751 756]); (2) x = 43; xlim([248 253]); (3) x = 38; xlim([1558 1563]);

beh = modAChDA; beh(42) = [];
fig = figure; fig.Position(3) = 1375;

sp(1) = subplot(1,3,1); hold on; 
for x = 43
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

plot(beh(x).time, fp_ach, 'g'); 
plot(beh(x).time, fp_da, 'm');
if ~isempty(beh(x).reward)
    stem(beh(x).reward./Fs, 10.*ones(length(beh(x).reward),1), 'b');
    lick = beh(x).lick(:)./Fs;
    a = [diff(lick.*1000) > 50]; % Identify licks that are <50ms after previous lick
    lick_sub = lick; lick_sub(1) = [];
    lick = [lick(1); lick_sub(a)];
    plot([lick, lick]',[10;11].*ones(2,length(lick)),'-k');
end
ylabel('Photometry (%dF/F)'); % ylim([-20 20])
acc = getAcc(beh(x).vel); 
plot(beh(x).time, acc - 5, 'k');
title(sprintf('%s',strtok(beh(x).rec,'-'))); axis square;

end
% xlim([80 85]);
xlim([751 756]);

sp(2) = subplot(1,3,2); hold on; 
for x = 43
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

plot(beh(x).time, fp_ach, 'g'); 
plot(beh(x).time, fp_da, 'm');
if ~isempty(beh(x).reward)
    stem(beh(x).reward./Fs, 10.*ones(length(beh(x).reward),1), 'b');
    lick = beh(x).lick(:)./Fs;
    a = [diff(lick.*1000) > 50]; % Identify licks that are <50ms after previous lick
    lick_sub = lick; lick_sub(1) = [];
    lick = [lick(1); lick_sub(a)];
    plot([lick, lick]',[10;11].*ones(2,length(lick)),'-k');
end
ylabel('Photometry (%dF/F)'); % ylim([-20 20])
acc = getAcc(beh(x).vel); 
plot(beh(x).time, acc - 5, 'k');
title(sprintf('%s',strtok(beh(x).rec,'-'))); axis square;

end
xlim([248 253]);

sp(3) = subplot(1,3,3); hold on; 
for x = 38
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

plot(beh(x).time, fp_ach, 'g'); 
plot(beh(x).time, fp_da, 'm');
if ~isempty(beh(x).reward)
    stem(beh(x).reward./Fs, 10.*ones(length(beh(x).reward),1), 'b');
    lick = beh(x).lick(:)./Fs;
    a = [diff(lick.*1000) > 50]; % Identify licks that are <50ms after previous lick
    lick_sub = lick; lick_sub(1) = [];
    lick = [lick(1); lick_sub(a)];
    plot([lick, lick]',[10;11].*ones(2,length(lick)),'-k');
end
ylabel('Photometry (%dF/F)'); % ylim([-20 20])
acc = getAcc(beh(x).vel); 
plot(beh(x).time, acc - 5, 'k');
title(sprintf('%s',strtok(beh(x).rec,'-'))); axis square;

end
xlim([1558 1563]);

linkaxes(sp,'y');

%% ACh DA - S1 (i) JM005
% 5 sec segments from (1) reward (2) locomotion (3) immobility
% 220506: (1) x = 43; xlim([751 756]); (2) x = 43; xlim([248 253]); (3) x = 38; xlim([1558 1563]);

beh = modAChDA; beh(42) = [];
fig = figure; fig.Position(3) = 1375;

sp(1) = subplot(1,3,1); hold on; 
for x = 34
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

plot(beh(x).time, fp_ach, 'g'); 
plot(beh(x).time, fp_da, 'm');
if ~isempty(beh(x).reward)
    stem(beh(x).reward./Fs, 10.*ones(length(beh(x).reward),1), 'b');
    lick = beh(x).lick(:)./Fs;
    a = [diff(lick.*1000) > 50]; % Identify licks that are <50ms after previous lick
    lick_sub = lick; lick_sub(1) = [];
    lick = [lick(1); lick_sub(a)];
    plot([lick, lick]',[10;11].*ones(2,length(lick)),'-k');
end
ylabel('Photometry (%dF/F)'); % ylim([-20 20])
acc = getAcc(beh(x).vel); 
plot(beh(x).time, acc - 5, 'k');
title(sprintf('%s',strtok(beh(x).rec,'-'))); axis square;

end
xlim([280 285]); % xlim([269 273]);

sp(2) = subplot(1,3,2); hold on; 
for x = 34
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

plot(beh(x).time, fp_ach, 'g'); 
plot(beh(x).time, fp_da, 'm');
if ~isempty(beh(x).reward)
    stem(beh(x).reward./Fs, 10.*ones(length(beh(x).reward),1), 'b');
    lick = beh(x).lick(:)./Fs;
    a = [diff(lick.*1000) > 50]; % Identify licks that are <50ms after previous lick
    lick_sub = lick; lick_sub(1) = [];
    lick = [lick(1); lick_sub(a)];
    plot([lick, lick]',[10;11].*ones(2,length(lick)),'-k');
end
ylabel('Photometry (%dF/F)'); % ylim([-20 20])
acc = getAcc(beh(x).vel); 
plot(beh(x).time, acc - 5, 'k');
title(sprintf('%s',strtok(beh(x).rec,'-'))); axis square;

end
xlim([1594 1599]);

sp(3) = subplot(1,3,3); hold on; 
for x = 34
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

plot(beh(x).time, fp_ach, 'g'); 
plot(beh(x).time, fp_da, 'm');
if ~isempty(beh(x).reward)
    stem(beh(x).reward./Fs, 10.*ones(length(beh(x).reward),1), 'b');
    lick = beh(x).lick(:)./Fs;
    a = [diff(lick.*1000) > 50]; % Identify licks that are <50ms after previous lick
    lick_sub = lick; lick_sub(1) = [];
    lick = [lick(1); lick_sub(a)];
    plot([lick, lick]',[10;11].*ones(2,length(lick)),'-k');
end
ylabel('Photometry (%dF/F)'); % ylim([-20 20])
acc = getAcc(beh(x).vel); 
plot(beh(x).time, acc - 5, 'k');
title(sprintf('%s',strtok(beh(x).rec,'-'))); axis square;

end
xlim([1671 1676]);

linkaxes(sp,'y');

%% ACh DA - S1 (ii) AK190
% 5 sec segments from (1) reward (2) locomotion (3) immobility
% 220506: (1) x = 43; xlim([751 756]); (2) x = 43; xlim([248 253]); (3) x = 38; xlim([1558 1563]);

beh = modAChDA; beh(42) = [];
fig = figure; fig.Position(3) = 1375;

sp(1) = subplot(1,3,1); hold on; 
for x = 44
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

plot(beh(x).time, fp_ach, 'g'); 
plot(beh(x).time, fp_da, 'm');
if ~isempty(beh(x).reward)
    stem(beh(x).reward./Fs, 10.*ones(length(beh(x).reward),1), 'b');
    lick = beh(x).lick(:)./Fs;
    a = [diff(lick.*1000) > 50]; % Identify licks that are <50ms after previous lick
    lick_sub = lick; lick_sub(1) = [];
    lick = [lick(1); lick_sub(a)];
    plot([lick, lick]',[10;11].*ones(2,length(lick)),'-k');
end
ylabel('Photometry (%dF/F)'); % ylim([-20 20])
acc = getAcc(beh(x).vel); 
plot(beh(x).time, acc - 5, 'k');
title(sprintf('%s',strtok(beh(x).rec,'-'))); axis square;

end
xlim([938 943]); % xlim([269 273]);

sp(2) = subplot(1,3,2); hold on; 
for x = 39
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

plot(beh(x).time, fp_ach, 'g'); 
plot(beh(x).time, fp_da, 'm');
if ~isempty(beh(x).reward)
    stem(beh(x).reward./Fs, 10.*ones(length(beh(x).reward),1), 'b');
    lick = beh(x).lick(:)./Fs;
    a = [diff(lick.*1000) > 50]; % Identify licks that are <50ms after previous lick
    lick_sub = lick; lick_sub(1) = [];
    lick = [lick(1); lick_sub(a)];
    plot([lick, lick]',[10;11].*ones(2,length(lick)),'-k');
end
ylabel('Photometry (%dF/F)'); % ylim([-20 20])
acc = getAcc(beh(x).vel); 
plot(beh(x).time, acc - 5, 'k');
title(sprintf('%s',strtok(beh(x).rec,'-'))); axis square;

end
xlim([22 27]);

sp(3) = subplot(1,3,3); hold on; 
for x = 39
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

plot(beh(x).time, fp_ach, 'g'); 
plot(beh(x).time, fp_da, 'm');
if ~isempty(beh(x).reward)
    stem(beh(x).reward./Fs, 10.*ones(length(beh(x).reward),1), 'b');
    lick = beh(x).lick(:)./Fs;
    a = [diff(lick.*1000) > 50]; % Identify licks that are <50ms after previous lick
    lick_sub = lick; lick_sub(1) = [];
    lick = [lick(1); lick_sub(a)];
    plot([lick, lick]',[10;11].*ones(2,length(lick)),'-k');
end
ylabel('Photometry (%dF/F)'); % ylim([-20 20])
acc = getAcc(beh(x).vel); 
plot(beh(x).time, acc - 5, 'k');
title(sprintf('%s',strtok(beh(x).rec,'-'))); axis square;

end
xlim([1091.5 1096.5]);

linkaxes(sp,'y');

%% ACh DA - S1 (iii) AK190
% 5 sec segments from (1) reward (2) locomotion (3) immobility
% 220506: (1) x = 43; xlim([751 756]); (2) x = 43; xlim([248 253]); (3) x = 38; xlim([1558 1563]);

beh = modAChDA; beh(42) = [];
fig = figure; fig.Position(3) = 1375;

sp(1) = subplot(1,3,1); hold on; 
for x = 41
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

plot(beh(x).time, fp_ach, 'g'); 
plot(beh(x).time, fp_da, 'm');
if ~isempty(beh(x).reward)
    stem(beh(x).reward./Fs, 10.*ones(length(beh(x).reward),1), 'b');
    lick = beh(x).lick(:)./Fs;
    a = [diff(lick.*1000) > 50]; % Identify licks that are <50ms after previous lick
    lick_sub = lick; lick_sub(1) = [];
    lick = [lick(1); lick_sub(a)];
    plot([lick, lick]',[10;11].*ones(2,length(lick)),'-k');
end
ylabel('Photometry (%dF/F)'); % ylim([-20 20])
acc = getAcc(beh(x).vel); 
plot(beh(x).time, acc - 5, 'k');
title(sprintf('%s',strtok(beh(x).rec,'-'))); axis square;

end
xlim([1598 1603]); % xlim([269 273]);

sp(2) = subplot(1,3,2); hold on; 
for x = 41
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

plot(beh(x).time, fp_ach, 'g'); 
plot(beh(x).time, fp_da, 'm');
if ~isempty(beh(x).reward)
    stem(beh(x).reward./Fs, 10.*ones(length(beh(x).reward),1), 'b');
    lick = beh(x).lick(:)./Fs;
    a = [diff(lick.*1000) > 50]; % Identify licks that are <50ms after previous lick
    lick_sub = lick; lick_sub(1) = [];
    lick = [lick(1); lick_sub(a)];
    plot([lick, lick]',[10;11].*ones(2,length(lick)),'-k');
end
ylabel('Photometry (%dF/F)'); % ylim([-20 20])
acc = getAcc(beh(x).vel); 
plot(beh(x).time, acc - 5, 'k');
title(sprintf('%s',strtok(beh(x).rec,'-'))); axis square;

end
xlim([329.5 334.5]);

sp(3) = subplot(1,3,3); hold on; 
for x = 41
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

plot(beh(x).time, fp_ach, 'g'); 
plot(beh(x).time, fp_da, 'm');
if ~isempty(beh(x).reward)
    stem(beh(x).reward./Fs, 10.*ones(length(beh(x).reward),1), 'b');
    lick = beh(x).lick(:)./Fs;
    a = [diff(lick.*1000) > 50]; % Identify licks that are <50ms after previous lick
    lick_sub = lick; lick_sub(1) = [];
    lick = [lick(1); lick_sub(a)];
    plot([lick, lick]',[10;11].*ones(2,length(lick)),'-k');
end
ylabel('Photometry (%dF/F)'); % ylim([-20 20])
acc = getAcc(beh(x).vel); 
plot(beh(x).time, acc - 5, 'k');
title(sprintf('%s',strtok(beh(x).rec,'-'))); axis square;

end
xlim([1581.5 1587.5]);

linkaxes(sp,'y');

%% ACh DA - S1 (iv) AK189
% 5 sec segments from (1) reward (2) locomotion (3) immobility
% 220506: (1) x = 43; xlim([751 756]); (2) x = 43; xlim([248 253]); (3) x = 38; xlim([1558 1563]);

beh = modAChDA; beh(42) = [];
fig = figure; fig.Position(3) = 1375;

sp(1) = subplot(1,3,1); hold on; 
for x = 41
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

plot(beh(x).time, fp_ach, 'g'); 
plot(beh(x).time, fp_da, 'm');
if ~isempty(beh(x).reward)
    stem(beh(x).reward./Fs, 10.*ones(length(beh(x).reward),1), 'b');
    lick = beh(x).lick(:)./Fs;
    a = [diff(lick.*1000) > 50]; % Identify licks that are <50ms after previous lick
    lick_sub = lick; lick_sub(1) = [];
    lick = [lick(1); lick_sub(a)];
    plot([lick, lick]',[10;11].*ones(2,length(lick)),'-k');
end
ylabel('Photometry (%dF/F)'); % ylim([-20 20])
acc = getAcc(beh(x).vel); 
plot(beh(x).time, acc - 5, 'k');
title(sprintf('%s',strtok(beh(x).rec,'-'))); axis square;

end
xlim([1598 1603]); % xlim([269 273]);

sp(2) = subplot(1,3,2); hold on; 
for x = 41
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

plot(beh(x).time, fp_ach, 'g'); 
plot(beh(x).time, fp_da, 'm');
if ~isempty(beh(x).reward)
    stem(beh(x).reward./Fs, 10.*ones(length(beh(x).reward),1), 'b');
    lick = beh(x).lick(:)./Fs;
    a = [diff(lick.*1000) > 50]; % Identify licks that are <50ms after previous lick
    lick_sub = lick; lick_sub(1) = [];
    lick = [lick(1); lick_sub(a)];
    plot([lick, lick]',[10;11].*ones(2,length(lick)),'-k');
end
ylabel('Photometry (%dF/F)'); % ylim([-20 20])
acc = getAcc(beh(x).vel); 
plot(beh(x).time, acc - 5, 'k');
title(sprintf('%s',strtok(beh(x).rec,'-'))); axis square;

end
xlim([329.5 334.5]);

sp(3) = subplot(1,3,3); hold on; 
for x = 41
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

plot(beh(x).time, fp_ach, 'g'); 
plot(beh(x).time, fp_da, 'm');
if ~isempty(beh(x).reward)
    stem(beh(x).reward./Fs, 10.*ones(length(beh(x).reward),1), 'b');
    lick = beh(x).lick(:)./Fs;
    a = [diff(lick.*1000) > 50]; % Identify licks that are <50ms after previous lick
    lick_sub = lick; lick_sub(1) = [];
    lick = [lick(1); lick_sub(a)];
    plot([lick, lick]',[10;11].*ones(2,length(lick)),'-k');
end
ylabel('Photometry (%dF/F)'); % ylim([-20 20])
acc = getAcc(beh(x).vel); 
plot(beh(x).time, acc - 5, 'k');
title(sprintf('%s',strtok(beh(x).rec,'-'))); axis square;

end
xlim([1581.5 1587.5]);

linkaxes(sp,'y');