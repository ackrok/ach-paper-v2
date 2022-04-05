%% Behavior
% propImm = []; % proportion of time spent immobile
% avgBoutDur = [];
% maxVel = [];
% avgVel = [];
% 
% for x = 1:length(beh)
%     imm = extractEventST([1:length(beh(x).vel)]',beh(x).onRest,beh(x).offRest,1);
%     propImm(x) = length(imm)/length(beh(x).vel);
%     Fs = 50;
%     avgBoutDur(x) = nanmean(beh(x).off - beh(x).on)/Fs;
%     mov = extractEventST([1:length(beh(x).vel)]',beh(x).on,beh(x).off,1);
%     maxVel(x) = max(beh(x).vel(mov));
%     tmp = [];
%     for y = 1:length(beh(x).on)
%         tmp(y) = nanmean(beh(x).vel(beh(x).on(y):beh(x).off(y)));
%     end
%     avgVel(x) = nanmean(tmp);
% end
load('C:\Users\Anya\Desktop\IV_LOCAL\data_comb\fullBeh_WT+6ohda_03-22.mat');

%%
figure; 
for x = 1:18
    subplot(6,3,x);
    plot(beh(x).vel);
end

%% compare WT vs lesion
fig = figure; fig.Position(3) = 1375;
subplot(1,3,1);
a = {[beh_wt.boutProp]*100, [beh_6ohda.boutProp]*100};
violinplot(a); xticklabels({'control','lesion'});
ylabel('Proportion of Time in Locomotion (%)'); ylim([0 100]);
title(sprintf('p = %1.3f',ranksum(a{1},a{2}))); axis('square');

subplot(1,3,2);
a = {[beh_wt.avgVel], [beh_6ohda.avgVel]};
violinplot(a); xticklabels({'control','lesion'});
ylabel('Average Velocity in Bout (cm/s)'); ylim([0 30]);
title(sprintf('p = %1.3f',ranksum(a{1},a{2}))); axis('square');

subplot(1,3,3);
a = {[beh_wt.avgBoutDur], [beh_6ohda.avgBoutDur]};
violinplot(a); 
xticklabels({'control','lesion'});
ylabel('Average Bout Duration (s)'); ylim([0 50]);
title(sprintf('p = %1.3f',ranksum(a{1},a{2}))); axis('square');
