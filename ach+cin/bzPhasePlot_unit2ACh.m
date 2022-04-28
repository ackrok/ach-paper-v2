[fName, fPath] = uigetfile('*.mat','MultiSelect','On');

%%
% beh = beh_da;
beh = behACh(1:13); % ACh recordings
grp = 1; % group labels 1-MSN, 2-CIN, 3-other

mat = struct;
for x = 1:length(beh)

    sub = WT(find([WT.label] == grp));
    ii = find(strcmp({sub.rec},beh(x).rec));
    if isempty(ii); continue; end
    
    if all(logical(~rem(beh(x).onRest,1))); d = 50; else; d = 1; end
    fp = struct; 
    fp.rec = beh(x).rec;
    fp.data = beh(x).FP{1}; % Extract photometry signal
    fp.channels = 1; fp.samplingRate = 50;
    
    st = cell(length(ii),1);
    for y = 1:length(ii)
        st{y} = extractEventST(sub(ii(y)).st, beh(x).onRest./d, beh(x).offRest./d, 1);
        z = find(st{y}.*fp.samplingRate > length(fp.data),1);
        if ~isempty(z); st{y} = st{y}(1:idx-1); end
    end
    spikes = struct; spikes.times = st;

    data = bz_PhaseModulation('spikes',spikes, 'lfp',fp, 'passband', [0.1 10], ...
        'intervals',[0 inf], 'samplingRate',50, 'method','hilbert', ...
        'numBins',30, 'powerThresh',0.5, 'plotting',false, 'saveMat',false);
    % fprintf('Spikes - FP done.\n');
    
    mat(x).rec = beh(x).rec;
    mat(x).n = [sub(ii).n];
    mat(x).phasebins = data.phasebins;
    mat(x).phasedistros = data.phasedistros;
    mat(x).pref = rad2deg(data.phasestats.m);
    mat(x).p = data.phasestats.p;
end
phasebins = data.phasebins;
%%
fig = figure; fig.Position(3) = 1375;
sm = 5;

switch grp
    case 2
        subplot(1,3,1); hold on % Example CIN
        x = find(strcmp({mat.rec},'IV066_rec01')); y = 4; phasedistros = movmean(mat(x).phasedistros(:,y),sm);
        plot([rad2deg(phasebins) 360+rad2deg(phasebins)], [phasedistros phasedistros].*100, 'b');
        x = find(strcmp({mat.rec},'IV043_rec02')); y = 1; phasedistros = movmean(mat(x).phasedistros(:,y),sm);
        plot([rad2deg(phasebins) 360+rad2deg(phasebins)], [phasedistros phasedistros].*100, 'c');
        x = find(strcmp({mat.rec},'IV066_rec03')); y = 3; phasedistros = movmean(mat(x).phasedistros(:,y),sm);
        plot([rad2deg(phasebins) 360+rad2deg(phasebins)], [phasedistros phasedistros].*100, 'g');
        ylim([1 6]);
    case 1
        subplot(1,3,1); hold on % Example MSN
        x = find(strcmp({mat.rec},'IV046_rec02')); y = 1; phasedistros = movmean(mat(x).phasedistros(:,y),sm);
        plot([rad2deg(phasebins) 360+rad2deg(phasebins)], [phasedistros phasedistros].*100, 'b');
        x = find(strcmp({mat.rec},'IV066_rec01')); y = 25; phasedistros = movmean(mat(x).phasedistros(:,y),sm);
        plot([rad2deg(phasebins) 360+rad2deg(phasebins)], [phasedistros phasedistros].*100, 'c');
        x = find(strcmp({mat.rec},'IV043_rec01')); y = 33; phasedistros = movmean(mat(x).phasedistros(:,y),sm);
        plot([rad2deg(phasebins) 360+rad2deg(phasebins)], [phasedistros phasedistros].*100, 'g');
        ylim([0 10]);
    case 3
        subplot(1,3,1); hold on % Example FSI
        x = find(strcmp({mat.rec},'IV066_rec01')); y = 1; phasedistros = movmean(mat(x).phasedistros(:,y),sm);
        plot([rad2deg(phasebins) 360+rad2deg(phasebins)], [phasedistros phasedistros].*100, 'b');
        x = find(strcmp({mat.rec},'IV066_rec03')); y = 5; phasedistros = movmean(mat(x).phasedistros(:,y),sm);
        plot([rad2deg(phasebins) 360+rad2deg(phasebins)], [phasedistros phasedistros].*100, 'c');
        x = find(strcmp({mat.rec},'IV043_rec01')); y = 2; phasedistros = movmean(mat(x).phasedistros(:,y),sm);
        plot([rad2deg(phasebins) 360+rad2deg(phasebins)], [phasedistros phasedistros].*100, 'g');
        ylim([1 6]);
end
yyaxis right
plot([rad2deg(phasebins) 360+rad2deg(phasebins)],  [cos(phasebins) cos(phasebins)], 'k');    
% title(sprintf('phase = %d, p = %1.3f',round(mat(x).pref(y)),mat(x).p(y)));
switch grp; case 1; lbl = 'pMSN'; case 2; lbl = 'pCIN'; case 3; lbl = 'pFSI'; end
title(sprintf('%s tuned to ACh photometry',lbl));
xlabel('phase'); xlim([0 720]); xticks([0:180:720]);
axis('square');

subplot(1,3,2); % Average across all units
phasedistros = movmean([mat.phasedistros],sm);
phasedistros(:,[mat.p] >= 0.05) = []; % remove untuned units
shadederrbar(rad2deg(phasebins), nanmean(phasedistros,2).*100, SEM(phasedistros,2).*100, 'b');
ylim([2 5]);
yyaxis right
plot(rad2deg(phasebins),  cos(phasebins), 'k');    
title(sprintf('n = %d units / %d are tuned',size(phasedistros,2),size([mat.phasedistros],2)));
xlabel('phase'); xlim([0 360]); xticks([0:90:360]);
axis('square');

subplot(1,3,3); % Phase preferences for all units
pref = [mat.pref]; pref([mat.p] >= 0.05) = []; 
rose(deg2rad(pref));
movegui(gcf,'center');