%%
fig = figure; fig.Position(3) = 1375;

subplot(1,3,1);
beh = modAChDA; x = 39;
plot(beh(x).time, beh(x).FP{2}); ylim([-5 15])
yyaxis right; plot(beh(x).time, getAcc(beh(x).vel),'k'); ylim([-1 1]);
title(sprintf('control %s',beh(x).rec)); axis('square');
xlim([1199.5 1209.5]); ylim([-5 15])

subplot(1,3,2);
beh = d2ant; x = 14;
plot(beh(x).time, beh(x).FP{1}); ylim([-5 15])
yyaxis right; plot(beh(x).time, getAcc(beh(x).vel),'k'); ylim([-1 1]);
title(sprintf('D2R antag %s',beh(x).rec)); axis('square');
xlim([410 420]);

subplot(1,3,3);
beh = behLes; x = 4;
plot(beh(x).time, beh(x).FP{2}); ylim([-5 15])
yyaxis right; plot(beh(x).time, getAcc(beh(x).vel),'k'); ylim([-1 1]);
title(sprintf('%s',beh(x).rec)); axis('square');
xlim([200 210]);

%%
figure;
beh = modAChDA; x = 39;
plot(beh(x).time, beh(x).FP{2});
yyaxis right; plot(beh(x).time, getAcc(beh(x).vel),'k'); ylim([-1 1]);
title(sprintf('%s',beh(x).rec));
xlim([1199.5 1209.5]); ylim([-5 15])

%%
% load('C:\Users\Anya\Desktop\FP_LOCAL\JM_DA_cw+D2Rant.mat')
% figure;
% beh = d2ant; 
% for x = 1:length(beh);
% subplot(4,4,x); hold on
% plot(beh(x).time, beh(x).FP{1});
% plot(beh(x).time, beh(x).FP{2});
% yyaxis right; plot(beh(x).time, getAcc(beh(x).vel),'k'); ylim([-1 1]);
% title(sprintf('%s',beh(x).rec));
% end

figure;
beh = d2ant; x = 14;
plot(beh(x).time, beh(x).FP{1}); ylim([-5 15])
yyaxis right; plot(beh(x).time, getAcc(beh(x).vel),'k'); ylim([-1 1]);
title(sprintf('%s',beh(x).rec));
xlim([410 420]);

%%
% load('R:\homes\ack466\ACh paper\7_6ohda\6ohda_fpdata.mat')

% figure;
% beh = behLes; 
% for x = 1:length(beh)
% subplot(4,3,x); hold on
% plot(beh(x).time, beh(x).FP{2});
% yyaxis right; plot(beh(x).time, getAcc(beh(x).vel),'k'); ylim([-1 1]);
% title(sprintf('%s',beh(x).rec));
% end

figure;
beh = behLes; x = 4;
plot(beh(x).time, beh(x).FP{2}); ylim([-5 15])
yyaxis right; plot(beh(x).time, getAcc(beh(x).vel),'k'); ylim([-1 1]);
title(sprintf('%s',beh(x).rec));
xlim([200 210]);