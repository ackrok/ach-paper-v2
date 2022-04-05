%%
load('C:\Users\Anya\Desktop\IV_LOCAL\data_comb\fullcin_6-2021_WT');

%% Firing rate: IMM vs MOV
sub = cinLes; beh = behLesion;
fr = [sub.fr];
fr_mov = []; fr_imm = []; 
for x = 1:length(sub)
    b = find(strcmp({beh.rec},sub(x).rec)); % Matching beh recording
    diffFs = 50; % diffFs = 50 if onset/offset are in samples
    fr_mov(x) = 1/mean(diff(extractEventST(sub(x).st, beh(b).on./diffFs, beh(b).off./diffFs, 0)));
    fr_imm(x) = 1/mean(diff(extractEventST(sub(x).st, beh(b).onRest./diffFs, beh(b).offRest./diffFs, 0)));
end

fr_wt = [];
for x = 1:length(cinWT)
    b = find(strcmp({behWT.rec},cinWT(x).rec)); % Matching beh recording
    diffFs = 1; % diffFs = 50 if onset/offset are in samples
    fr_wt(2,x) = 1/mean(diff(extractEventST(cinWT(x).st, behWT(b).on./diffFs, behWT(b).off./diffFs, 0)));
    fr_wt(1,x) = 1/mean(diff(extractEventST(cinWT(x).st, behWT(b).onRest./diffFs, behWT(b).offRest./diffFs, 0)));
end

%% PLOT firing rate: IMM vs MOV
figure; hold on
a = [fr_imm;fr_mov];
plot(a,'--k');
plot(nanmean(a,2),'.r','MarkerSize',20); errorbar(nanmean(a,2), SEM(a,2), 'r');
plot(nanmean(fr_wt,2),'.k','MarkerSize',20); errorbar(nanmean(fr_wt,2), SEM(a,2), 'k');
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'imm','mov'});
ylabel('Firing Rate (Hz)'); ylim([2 10]); yticks([0:2:10]);
title(sprintf('CIN FR lesion imm(%1.2f) mov(%1.2f) p = %1.3f \n WT imm(%1.2f) mov(%1.2f) p = %1.3f \n WT vs LES imm(p %1.3f) mov(p %1.3f)',...
    nanmean(fr_imm),nanmean(fr_mov),signrank(fr_imm,fr_mov),...
    nanmean(fr_wt(1,:)),nanmean(fr_wt(2,:)),signrank(fr_wt(1,:),fr_wt(2,:)),...
    ranksum(fr_imm,fr_wt(1,:)), ranksum(fr_mov,fr_wt(2,:))));
axis('square');

%% CIN activity to acceleration locs
