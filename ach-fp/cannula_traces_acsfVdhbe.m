%%
fig = figure; fig.Position(3) = 1375;
load('C:\Users\Anya\Desktop\FP_LOCAL\AK189-190_cannula.mat')
beh = cannula;
x = 21;

idx_imm = extractEventST([1:length(fp)]', beh(x).onRest, beh(x).offRest, 1); % Event times during rest
bin_imm = zeros(length(fp),1); bin_imm(idx_imm) = 1;

subplot(1,3,1); hold on
area(beh(x).time, 25.*bin_imm,'FaceAlpha',0.2,'FaceColor','r','EdgeAlpha',0);
plot(beh(x).time, beh(x).FP{1}, 'Color', [0.05 0.75 0.45]);
plot(beh(x).time, beh(x).FP{2}, 'm');
plot(beh(x).time, getAcc(beh(x).vel), 'k');
title(sprintf('%s',beh(x).rec))
xlim([1923 1933]);
axis('square');

%
load('C:\Users\Anya\Desktop\FP_LOCAL\AK189-197_cannula+rew_v4.mat')
x = 1; y = 4;
beh = s(y).s;

idx_imm = extractEventST([1:length(fp)]', beh(x).onRest, beh(x).offRest, 1); % Event times during rest
bin_imm = zeros(length(fp),1); bin_imm(idx_imm) = 1;

subplot(1,3,2); hold on
area(beh(x).time, 25.*bin_imm,'FaceAlpha',0.2,'FaceColor','r','EdgeAlpha',0);
plot(beh(x).time, beh(x).FP{1}, 'Color', [0.05 0.75 0.45]);
plot(beh(x).time, beh(x).FP{2}, 'm');
plot(beh(x).time, getAcc(beh(x).vel), 'k');
title(sprintf('%s',beh(x).rec))
xlim([1940 1950]); xticks([1:6000]);
axis('square');

%
load('C:\Users\Anya\Desktop\FP_LOCAL\AK189-197_cannula+rew_v4.mat')
x = 2; y = 4;
beh = s(y).s;

idx_imm = extractEventST([1:length(fp)]', beh(x).onRest, beh(x).offRest, 1); % Event times during rest
bin_imm = zeros(length(fp),1); bin_imm(idx_imm) = 1;
subplot(1,3,3); hold on
area(beh(x).time, 25.*bin_imm,'FaceAlpha',0.2,'FaceColor','r','EdgeAlpha',0);
plot(beh(x).time, beh(x).FP{1}, 'Color', [0.05 0.75 0.45]);
plot(beh(x).time, beh(x).FP{2}, 'm');
plot(beh(x).time, getAcc(beh(x).vel), 'k');
title(sprintf('%s',beh(x).rec))
xlim([3155 3165]);
axis('square');