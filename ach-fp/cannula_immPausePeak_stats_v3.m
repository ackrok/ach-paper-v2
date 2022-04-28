z = 2; % PAUSE(1) or PEAK(2)

switch z; case 1; lbl = 'Pause'; case 2; lbl = 'Peak'; end
nAn = 4;
%% PLOT
fig = figure; fig.Position(3) = 1375;

%IMM - frequency
freq_mat = [];
for y = 1:3; freq_mat(:,y) = s(y).freq(:,z); end
freq_mat(isnan(freq_mat)) = 0;
a = freq_mat./nanmean(freq_mat(:,1));

subplot(1,3,1); hold on
plot([1 3],[1 1],'--k');
errorbar(abs(nanmean(a)'),abs(SEM(a,1)'),'.k','MarkerSize',20);
plot([1.15; 1.85; 2.85].*ones(3,nAn),abs(a'),'.k','MarkerSize',20);
plot([1.15; 1.85].*ones(2,nAn),abs(a(:,[1 2])'),':b','MarkerSize',20);
plot([1.15; 2.85].*ones(2,nAn),abs(a(:,[1 3])'),':r','MarkerSize',20);
xlim([0.5 3.5]); xticks([1:3]); xticklabels({'aCSF','D1R/D2R','NMDA/AMPA'});
ylabel(sprintf('%s Frequency (Hz)',lbl)); ylim([0 2]); yticks([0:0.5:3]);
p = []; for x = 2:3; [~, p(x)] = ttest(a(:,1),a(:,x)); end
title(sprintf('IMM: %s Freq (a/da: %1.3f)(a/glu: %1.3f)',lbl,p(2),p(3))); axis('square');

%IMM - duration
dur_mat = [];
for y = 1:3; for x = 1:4; dur_mat(x,y) = nanmean(s(y).dur{x,z}); end; end
dur_mat(isnan(dur_mat)) = 0;
dur_mat = dur_mat.*(1000/Fs); % Adjust from samples to ms
r = [1 2]; a = dur_mat(:,r)./nanmean(dur_mat(:,1));

subplot(1,3,2); hold on
plot([1 2],[1 1],'--k');
errorbar(abs(nanmean(a)'),abs(SEM(a,1)'),'.k','MarkerSize',20);
plot([1.15; 1.85].*ones(2,nAn),abs(a'),'.k','MarkerSize',20);
plot([1.15; 1.85].*ones(2,nAn),abs(a(:,[1 2])'),':b','MarkerSize',20);
% plot([1.15; 2.85].*ones(2,nAn),abs(a(:,[1 3])'),':r','MarkerSize',20);
xlim([0.5 2.5]); xticks([1:2]); xticklabels({'aCSF','D1R/D2R'});
ylabel(sprintf('%s Duration',lbl)); ylim([0 2]); yticks([0:0.5:3]);
[~, p] = ttest(a(:,1),a(:,2));
title(sprintf('IMM: %s Dur (a/da: %1.3f)',lbl,p)); axis('square');

%IMM - amplitude
amp_mat = [];
for y = 1:3; for x = 1:4; amp_mat(x,y) = nanmean(s(y).amp{x,z}); end; end
amp_mat(isnan(amp_mat)) = 0;
switch z; case 1; amp_mat = -1*amp_mat; end
r = [1 2]; a = amp_mat(:,r)./nanmean(amp_mat(:,1));

subplot(1,3,3); hold on
plot([1 2],[1 1],'--k');
errorbar(abs(nanmean(a)'),abs(SEM(a,1)'),'.k','MarkerSize',20);
plot([1.15; 1.85].*ones(2,nAn),abs(a'),'.k','MarkerSize',20);
plot([1.15; 1.85].*ones(2,nAn),abs(a(:,[1 2])'),':b','MarkerSize',20);
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'aCSF','D1R/D2R'});
ylabel(sprintf('%s Amplitude',lbl)); ylim([0 2]); yticks([0:0.5:3]);
[~, p] = ttest(a(:,1),a(:,2));
title(sprintf('IMM: %s Amp (a/da: %1.3f)',lbl,p)); axis('square');

movegui(gcf,'center');