z = 2; % PAUSE(1) or PEAK(2)

switch z; case 1; lbl = 'Pause'; case 2; lbl = 'Peak'; end
nAn = 8;
%% PLOT
fig = figure; fig.Position(3) = 1375;

%IMM - frequency
freq_mat = s(1).freq(:,z); %[s(1).freq([1:4],z); s(1).freq([1:4],z)];
freq_mat(:,2) = [s(2).freq([1:4],z); nan; s(2).freq(5,z); nan; s(2).freq(6,z)];
freq_mat(:,3) = s(3).freq(:,z);
switch z; case 1; freq_mat([1:3,6,7],3) = 0; case 2; freq_mat([3,4,6,7],3) = 0; end
freq_mat(:,4) = [s(2).freq([1:2],z); nan; s(2).freq(3,z); nan; s(2).freq(4,z); nan; s(2).freq(5,z)];
a = freq_mat./nanmean(freq_mat(:,1));
%
subplot(1,3,1); hold on
plot([1 4],[1 1],'--k');
errorbar(abs(nanmean(a)'),abs(SEM(a,1)'),'.k','MarkerSize',20);
plot([1.15; 1.85; 2.85; 3.85].*ones(4,nAn),abs(a'),'.k','MarkerSize',20);
plot([1.15; 1.85].*ones(2,nAn),abs(a(:,[1 2])'),':m','MarkerSize',20);
plot([1.15; 2.85].*ones(2,nAn),abs(a(:,[1 3])'),':r','MarkerSize',20);
plot([1.15; 3.85].*ones(2,nAn),abs(a(:,[1 4])'),':b','MarkerSize',20);
xlim([0.5 4.5]); xticks([1:4]); xticklabels({'aCSF','D1R/D2R','NMDA/AMPA','DHbE'});
ylabel(sprintf('IMM %s Frequency',lbl)); ylim([0 3]); yticks([0:0.5:3]);
% p = []; for x = 2:3; [~, p(x)] = ttest(a(:,1),a(:,x)); end
% [~,p(4)] = ttest(a([1 2 4],1),a([1 2 4],4));
p = []; for x = 2:4; p(x) = signrank(a(:,1),a(:,x)); end
title(sprintf('(a/da: %1.4f)(a/glu: %1.4f)(a/n: %1.4f)',p(2),p(3),p(4))); axis('square');

%IMM - duration
dur_mat = [];
y = 1; tmp = cellfun(@nanmean, s(y).dur); dur_mat(:,y) = tmp(:,z);
y = 2; tmp = cellfun(@nanmean, s(y).dur); dur_mat(:,y) = [tmp([1:4],z); nan; tmp(5,z); nan; tmp(6,z)];
y = 3; dur_mat(:,3) = nan(nAn,1);
y = 4; tmp = cellfun(@nanmean, s(y).dur); dur_mat(:,y) = [tmp([1:2],z); nan; tmp(3,z); nan; tmp(4,z); nan; tmp(5,z)];
dur_mat = dur_mat.*(1000/Fs); % Adjust from samples to ms
a = dur_mat./nanmean(dur_mat(:,1));

subplot(1,3,2); hold on
plot([1 3],[1 1],'--k');
errorbar(abs(nanmean(a)'),abs(SEM(a,1)'),'.k','MarkerSize',20);
plot([1.15; 1.85; 2.85; 3.85].*ones(4,nAn),abs(a'),'.k','MarkerSize',20);
plot([1.15; 1.85].*ones(2,nAn),abs(a(:,[1 2])'),':m','MarkerSize',20);
plot([1.15; 3.85].*ones(2,nAn),abs(a(:,[1 4])'),':b','MarkerSize',20);
xlim([0.5 4.5]); xticks([1:3]); xticklabels({'aCSF','D1R/D2R','glu','DHbE'});
ylabel(sprintf('IMM %s Duration',lbl)); ylim([0 2.5]); yticks([0:0.5:3]);
% [~, p] = ttest(a(:,1),a(:,2)); [~,p(2)] = ttest(a([1 2 4],1),a([1 2 4],3));
p = []; for x = [2 4]; p(x) = signrank(a(:,1),a(:,x)); end
title(sprintf('(a/da: %1.4f)(a/n: %1.4f)',p(2),p(4))); axis('square');

%IMM - amplitude
amp_mat = [];
y = 1; tmp = cellfun(@nanmean, s(y).amp); amp_mat(:,y) = tmp(:,z);
y = 2; tmp = cellfun(@nanmean, s(y).amp); amp_mat(:,y) = [tmp([1:4],z); nan; tmp(5,z); nan; tmp(6,z)];
y = 3; amp_mat(:,3) = nan(nAn,1);
y = 4; tmp = cellfun(@nanmean, s(y).amp); amp_mat(:,y) = [tmp([1:2],z); nan; tmp(3,z); nan; tmp(4,z); nan; tmp(5,z)];
switch z; case 1; amp_mat = -1*amp_mat; end
a = amp_mat./nanmean(amp_mat(:,1));

subplot(1,3,3); hold on
plot([1 3],[1 1],'--k');
errorbar(abs(nanmean(a)'),abs(SEM(a,1)'),'.k','MarkerSize',20);
plot([1.15; 1.85; 2.85; 3.85].*ones(4,nAn),abs(a'),'.k','MarkerSize',20);
plot([1.15; 1.85].*ones(2,nAn),abs(a(:,[1 2])'),':m','MarkerSize',20);
plot([1.15; 3.85].*ones(2,nAn),abs(a(:,[1 4])'),':b','MarkerSize',20);
xlim([0.5 4.5]); xticks([1:3]); xticklabels({'aCSF','D1R/D2R','glu','DHbE'});
ylabel(sprintf('IMM %s Amplitude',lbl)); ylim([0 2.5]); yticks([0:0.5:3]);
% [~, p] = ttest(a(:,1),a(:,2)); [~,p(2)] = ttest(a([1 2 4],1),a([1 2 4],3));
p = []; for x = [2 4]; p(x) = signrank(a(:,1),a(:,x)); end
title(sprintf('(a/da: %1.4f)(a/n: %1.4f)',p(2),p(4))); axis('square');

movegui(gcf,'center');