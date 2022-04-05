%% UNIT 3D PLOT -- Peak Latency (ms), Prop (ISI > 2s), Firing Rate (Hz)
roi = lesion; roi2 = cinLes;
%
peakL = [roi.peakL]; 
fr_les    = [roi.fr]; 
pISI2 = [roi.pISI2];
%
clr = {'k','b','g'};
figure; hold on
for x = [1 3]
    idx = find([roi.label] == x);
    %scatter3(pISI2(idx), peakL(idx), fr(idx), 'DisplayName', sprintf('Type %d',x));
    plot3(pISI2(idx),peakL(idx),fr_les(idx),'.','MarkerSize',15,'Color',clr{x},'DisplayName', sprintf('Type %d',x));
end
x = 2; plot3([roi2.pISI2],[roi2.peakL],[roi2.fr],'.','MarkerSize',15,'Color',clr{x},'DisplayName', sprintf('Type %d',x));

xlabel('prop (ISI) > 2s'); ylabel('Peak Latency (ms)'); zlabel('Firing Rate (Hz)');
title(sprintf('Control Units (n = %d), t1-%d, t2-%d, t3-%d',length(roi),length(find([roi.label] == 1)),length(find([roi.label] == 2)),length(find([roi.label] == 3))));
legend; grid on; xlim([0 1]); ylim([0 3]); zlim([0 50]); view(45,45); 

%% Plot Average Waveforms
roi = lesion;

fig = figure; fig.Position(3) = 1375;
sp(1) = subplot(1,3,1); hold on
wf_time = [-4:1/30:4];
sub = lesion(find([lesion.label] == 1)); wf = [sub.wf]; wf = zscore(wf);
shadederrbar(wf_time,nanmean(wf,2),SEM(wf,2),'m');
sub = WT(find([WT.label] == 1)); wf = [sub.wf]; wf = zscore(wf);
plot(wf_time, nanmean(wf,2), 'k');
xlabel('Time (ms)'); ylabel('Voltage (uV)'); title('SPN'); axis('square');

sp(2) = subplot(1,3,2); hold on
sub = lesion(find([lesion.label] == 2)); wf = [sub.wf]; wf = zscore(wf);
a = [];
for x = 1:length(sub)
    [m,ii] = max(abs(sub(x).wf));
    if sub(x).wf(ii) < 0; a = [a,x]; end
end
shadederrbar(wf_time,nanmean(wf(:,a),2),SEM(wf(:,a),2),'b');
sub = WT(find([WT.label] == 2)); wf = [sub.wf]; wf = zscore(wf);
plot(wf_time, nanmean(wf,2), 'k');
xlabel('Time (ms)'); ylabel('Voltage (uV)'); title('CIN'); axis('square');

sp(3) = subplot(1,3,3); hold on
sub = lesion(find([lesion.label] == 3)); wf = [sub.wf]; wf = zscore(wf);
shadederrbar(wf_time,nanmean(wf,2),SEM(wf,2),'g');
sub = WT(find([WT.label] == 3)); wf = [sub.wf]; wf = zscore(wf);
plot(wf_time, nanmean(wf,2), 'k');
xlabel('Time (ms)'); ylabel('Voltage (uV)'); title('FSI'); axis('square');
linkaxes(sp,'y'); linkaxes(sp,'x'); xlim([-2 2])
movegui(gcf,'center');

%% CIN Firing Rate + ACG
st = {cinLes.st};

acgLin = []; %initialize matrixes
for x = 1:length(st)
    [xLin, nLin, ~, ~]= myACG(st{x},[],[]); %compute autocorrelogram with linear bins
    nLin = [flip(nLin);nLin]; %symmetrical ACG
    nLin = nLin./length(st{x});
    xLin = [-1.*flip(xLin);xLin]; xLin = xLin.*1000; %x-values for ACG in milliseconds
    acgLin = [acgLin, nLin]; %add to matrix
end

fig = figure; fig.Position(3) = 1375;
subplot(1,3,1); hold on
fr_les = [cinLes.fr];
[~, edges] = histcounts(log10(fr_les),'BinWidth',0.1); 
histogram(fr_les,10.^edges,'Normalization','probability'); %histogram of ISI distribution
fr_wt = [WT(find([WT.label] == 2)).fr];
[~, edges] = histcounts(log10(fr_wt),'BinWidth',0.1); 
histogram(fr_wt,10.^edges,'Normalization','probability'); %histogram of ISI distribution
set(gca,'xscale','log'); %set log-scale for x-axis
xlabel('Firing Rate (Hz)'); xlim([0.1 100]); % xticks([10^-1 10^0 10^2])
legend({'lesion','WT'});
title(sprintf('CIN FR - WT(%1.2f) LES(%1.2f) - p = %1.3f',nanmean(fr_wt),nanmean(fr_les),ranksum(fr_wt,fr_les))); axis('square');

subplot(1,3,2); hold on
cv_les = [cinLes.CV];
cv_wt = [WT(find([WT.label] == 2)).CV];
violinplot({cv_les,cv_wt}); xticklabels({'lesion','WT'});
ylabel('CV'); ylim([0 2.5]);
title(sprintf('CIN FR - WT(%1.2f) LES(%1.2f) - p = %1.3f',nanmean(cv_wt),nanmean(cv_les),ranksum(cv_les,cv_wt))); axis('square');

subplot(1,3,3);
sm = 10;
shadederrbar(xLin, movmean(nanmean(acgLin,2),sm), movmean(SEM(acgLin,2),sm), 'b');
ylabel('Firing Rate (Hz)');
xlabel('Lag (ms)'); xlim([-250 250]); 
title('CIN ACG'); axis('square');

%% MSN Firing Rate + ACG
st = {lesion(find([lesion.label] == 1)).st};

acgLin = []; %initialize matrixes
for x = 1:length(st)
    [xLin, nLin, ~, ~]= myACG(st{x},[],[]); %compute autocorrelogram with linear bins
    nLin = [flip(nLin);nLin]; %symmetrical ACG
    nLin = nLin./length(st{x});
    xLin = [-1.*flip(xLin);xLin]; xLin = xLin.*1000; %x-values for ACG in milliseconds
    acgLin = [acgLin, nLin]; %add to matrix
end

fig = figure; fig.Position(3) = 1000;
subplot(1,2,1); hold on
fr_les = [lesion(find([lesion.label] == 1)).fr];
[~, edges] = histcounts(log10(fr_les),'BinWidth',0.1); 
histogram(fr_les,10.^edges,'Normalization','probability'); %histogram of ISI distribution
fr_wt = [WT(find([WT.label] == 1)).fr];
[~, edges] = histcounts(log10(fr_wt),'BinWidth',0.1); 
histogram(fr_wt,10.^edges,'Normalization','probability'); %histogram of ISI distribution
set(gca,'xscale','log'); %set log-scale for x-axis
xlim([0.01 100]); % xticks([10^-1 10^0 10^2])
legend({'lesion','WT'});
title(sprintf('MSN FR - WT(%1.2f) LES(%1.2f) - p = %1.3f',nanmean(fr_wt),nanmean(fr_les),ranksum(fr_wt,fr_les))); axis('square');

subplot(1,2,2);
sm = 10;
shadederrbar(xLin, movmean(nanmean(acgLin,2),sm), movmean(SEM(acgLin,2),sm), 'b');
ylabel('Firing Rate (Hz)');
xlabel('Lag (ms)'); xlim([-250 250]); 
title('MSN ACG'); axis('square');
