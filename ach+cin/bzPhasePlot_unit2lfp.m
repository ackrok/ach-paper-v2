[fName, fPath] = uigetfile('*.mat','MultiSelect','On');

%%
sub = cinWT; beh = behWT;
% sub = WT(find([WT.label] == 1)); beh = behwt;

mat = struct;
h = waitbar(0, 'tuning unit to LFP');
for x = 1:length(fName)
    load(fullfile(fPath,fName{x}));

    ii = find(strcmp({sub.rec},lfp.rec));
    if isempty(ii); continue; end
    st = cell(length(ii),1);
    ii_b = find(strcmp({beh.rec},lfp.rec));
    if all(logical(~rem(beh(ii_b).onRest,1))); d = 50; else; d = 1; end
    for y = 1:length(ii)
        st{y} = extractEventST(sub(ii(y)).st, beh(ii_b).onRest./d, beh(ii_b).offRest./d, 1);
    end
    spikes = struct; spikes.times = st;

    data = bz_PhaseModulation('spikes',spikes, 'lfp',lfp, 'passband', [0.5 4], ...
        'intervals',[0 inf], 'samplingRate',1250, 'method','hilbert', ...
        'numBins',30, 'powerThresh',1, 'plotting',false, 'saveMat',false);
    % fprintf('Spikes - LFP done.\n');
    
    mat(x).rec = lfp.rec;
    mat(x).n = [sub(ii).n];
    mat(x).phasebins = data.phasebins;
    mat(x).phasedistros = data.phasedistros;
    mat(x).pref = rad2deg(data.phasestats.m);
    mat(x).p = data.phasestats.p;
    mat(x).kappa = data.phasestats.k;
    waitbar(x/length(fName),h);
end
phasebins = data.phasebins;
close(h);

%%
fig = figure; fig.Position(3) = 1375;
sm = 5;

% subplot(1,3,1); hold on % Example CIN
% x = find(strcmp({mat.rec},'IV066_rec01')); y = 2; 
% phasedistros = movmean(mat(x).phasedistros(:,y),sm);
% phasedistros = normalize(phasedistros,'range');
% plot([rad2deg(phasebins) 360+rad2deg(phasebins)], [phasedistros phasedistros], 'c');
% x = find(strcmp({mat.rec},'IV052_rec01')); y = 2; 
% phasedistros = movmean(mat(x).phasedistros(:,y),sm);
% phasedistros = normalize(phasedistros,'range');
% plot([rad2deg(phasebins) 360+rad2deg(phasebins)], [phasedistros phasedistros], 'b');
% x = find(strcmp({mat.rec},'IV046_rec02')); y = 1; 
% phasedistros = movmean(mat(x).phasedistros(:,y),sm);
% phasedistros = normalize(phasedistros,'range');
% plot([rad2deg(phasebins) 360+rad2deg(phasebins)], [phasedistros phasedistros], 'g');
% title('pCIN tuned to LFP');
% ylim([0 1]);

subplot(1,3,1); hold on % Example MSN
x = find(strcmp({mat.rec},'IV042_rec03')); y = 4; 
phasedistros = movmean(mat(x).phasedistros(:,y),sm);
phasedistros = normalize(phasedistros,'range');
plot([rad2deg(phasebins) 360+rad2deg(phasebins)], [phasedistros phasedistros], 'b');
x = find(strcmp({mat.rec},'IV066_rec01')); y = 25; 
phasedistros = movmean(mat(x).phasedistros(:,y),sm);
phasedistros = normalize(phasedistros,'range');
plot([rad2deg(phasebins) 360+rad2deg(phasebins)], [phasedistros phasedistros], 'c');
x = find(strcmp({mat.rec},'IV043_rec01')); y = 33; 
phasedistros = movmean(mat(x).phasedistros(:,y),sm);
phasedistros = normalize(phasedistros,'range');
plot([rad2deg(phasebins) 360+rad2deg(phasebins)], [phasedistros phasedistros], 'g');
ylim([0 1]);
title('pMSN tuned to LFP');

yyaxis right
plot([rad2deg(phasebins) 360+rad2deg(phasebins)],  [cos(phasebins) cos(phasebins)], 'k');    
% title(sprintf('phase = %d, p = %1.3f',round(mat(x).pref(y)),mat(x).p(y)));
xlabel('phase'); xlim([0 720]); xticks([0:180:720]);
axis('square');
%%
subplot(1,3,2); % Average across all units
phasedistros = movmean([mat.phasedistros],sm);
phasedistros(:,[mat.p] >= 0.05) = []; % remove untuned units
phasedistros = normalize(phasedistros,'range');
shadederrbar(rad2deg(phasebins), nanmean(phasedistros,2), SEM(phasedistros,2), 'b');
ylim([0 1]);
yyaxis right
plot(rad2deg(phasebins),  cos(phasebins), 'k');    
title(sprintf('n = %d units / %d are tuned',size(phasedistros,2),size([mat.phasedistros],2)));
xlabel('phase'); xlim([0 360]); xticks([0:90:360]);
axis('square');

subplot(1,3,3); % Phase preferences for all units
pref = [mat.pref]; pref([mat.p] >= 0.05) = []; 
rose(deg2rad(pref));
movegui(gcf,'center');