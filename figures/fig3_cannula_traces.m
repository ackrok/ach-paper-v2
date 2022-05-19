load('C:\Users\Anya\Desktop\FP_LOCAL\AK189-197_cannula+rew_full.mat')
beh = cannula; 

%% aCSF: AK190-210902
% load('C:\Users\Anya\Desktop\FP_LOCAL\AK189-190_cannula.mat')
beh = cannula;
x = 21;
figure;
sp(1) = subplot(2,1,1); hold on
plot(beh(x).time, beh(x).FP{1}, 'Color', [0.05 0.75 0.45]);
plot(beh(x).time, beh(x).FP{2}, 'm');
title(sprintf('%s - %s',beh(x).rec,beh(x).rx),'Interpreter','none');
ylabel('ACh (%dF/F)');ylim([-5 22])

sp(2) = subplot(2,1,2); hold on
    % idx_rest = extractEventST([1:length(beh(x).time)]',beh(x).onRest,beh(x).offRest,1);
    % bin_rest = zeros(length(beh(x).time),1); bin_rest(idx_rest) = 1;
    % area(beh(x).time, bin_rest, 'FaceColor','r', 'FaceAlpha',0.2);
    % idx_rew = extractEventST([1:length(beh(x).time)]',beh(x).reward,beh(x).reward+50,1);
    % bin_rew = zeros(length(beh(x).time),1); bin_rew(idx_rew) = 1;
    % area(beh(x).time, bin_rew, 'FaceColor','b', 'FaceAlpha',0.2);
plot(beh(x).time, getAcc(beh(x).vel), 'k');
ylabel('Acceleration (cm/s^2)');ylim([-1 1]);
linkaxes(sp,'x');
xlim([1922.5 1932.5]); 
xlim([1560 1570])
xlim([2185 2195])

%% iGluR: AK190-211104
x = 20;
figure;
sp(1) = subplot(2,1,1);
plot(beh(x).time, beh(x).FP{1}, 'Color', [0.85 0.35 0.1]);
title(sprintf('%s - %s',beh(x).rec,beh(x).rx),'Interpreter','none');
ylabel('ACh (%dF/F)'); ylim([-5 22])
 
sp(2) = subplot(2,1,2); hold on
    % idx_rest = extractEventST([1:length(beh(x).time)]',beh(x).onRest,beh(x).offRest,1);
    % bin_rest = zeros(length(beh(x).time),1); bin_rest(idx_rest) = 1;
    % area(beh(x).time, bin_rest, 'FaceColor','r', 'FaceAlpha',0.2);
    % idx_rew = extractEventST([1:length(beh(x).time)]',beh(x).reward,beh(x).reward+50,1);
    % bin_rew = zeros(length(beh(x).time),1); bin_rew(idx_rew) = 1;
    % area(beh(x).time, bin_rew, 'FaceColor','b', 'FaceAlpha',0.2);
plot(beh(x).time, getAcc(beh(x).vel), 'k');
ylabel('Acceleration (cm/s^2)'); ylim([-1 1]);
linkaxes(sp,'x');
xlim([3125 3135]);
 
%% d1d2: AK190-211110
x = 28; 
figure;
sp(1) = subplot(2,1,1);
plot(beh(x).time, beh(x).FP{1}, 'Color', [0.5 0.2 0.55]);
title(sprintf('%s - %s',beh(x).rec,beh(x).rx),'Interpreter','none');
ylabel('ACh (%dF/F)'); ylim([-5 22])
 
sp(2) = subplot(2,1,2); hold on
    idx_rest = extractEventST([1:length(beh(x).time)]',beh(x).onRest,beh(x).offRest,1);
    bin_rest = zeros(length(beh(x).time),1); bin_rest(idx_rest) = 1;
    area(beh(x).time, bin_rest, 'FaceColor','r', 'FaceAlpha',0.2);
    idx_rew = extractEventST([1:length(beh(x).time)]',beh(x).reward,beh(x).reward+50,1);
    bin_rew = zeros(length(beh(x).time),1); bin_rew(idx_rew) = 1;
    area(beh(x).time, bin_rew, 'FaceColor','b', 'FaceAlpha',0.2);
plot(beh(x).time, getAcc(beh(x).vel), 'k');
ylabel('Acceleration (cm/s^2)'); ylim([-1 1]);
linkaxes(sp,'x');
xlim([1190 1200]);

%% nAChR antagonist
% beh = s(4).s; x = 3;
figure;
sp(1) = subplot(2,1,1); hold on
plot(beh(x).time, beh(x).FP{1} - nanmean(beh(x).FP{1}), 'Color', [0.05 0.75 0.45]);
% plot(beh(x).time, beh(x).FP{2} - nanmean(beh(x).FP{2}), 'm');
title(sprintf('%s - %s',beh(x).rec,beh(x).rx),'Interpreter','none');
ylabel('ACh (%dF/F)'); ylim([-5 22])

sp(2) = subplot(2,1,2); hold on
    % idx_rest = extractEventST([1:length(beh(x).time)]',beh(x).onRest,beh(x).offRest,1);
    % bin_rest = zeros(length(beh(x).time),1); bin_rest(idx_rest) = 1;
    % area(beh(x).time, bin_rest, 'FaceColor','r', 'FaceAlpha',0.2);
    % idx_rew = extractEventST([1:length(beh(x).time)]',beh(x).reward,beh(x).reward+50,1);
    % bin_rew = zeros(length(beh(x).time),1); bin_rew(idx_rew) = 1;
    % area(beh(x).time, bin_rew, 'FaceColor','b', 'FaceAlpha',0.2);
plot(beh(x).time, getAcc(beh(x).vel), 'k');
ylabel('Acceleration (cm/s^2)'); ylim([-1 1]);
linkaxes(sp,'x');
xlim([3405 3415])
% xlim([2555 2565])