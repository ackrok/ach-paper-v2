%% Input Variables
load('C:\Users\Anya\Desktop\IV_LOCAL\data_comb\fullcin_2022-Feb_WT');
beh = behWT([11:22, 32:45]); % DLS AP orientation recordings only
% sub = cinWT([mat.n]);

%% Align acceleration to CIN unit spikes
Fs = 50; win = [-1 1]; %CHANGE: window for STA
align = {length(sub),2}; nShuff = 10; alignZ = align;
h = waitbar(0, 'STA: acceleration to CIN unit spike times');
for x = 1:length(sub)
    ib = find(strcmp({beh.rec},sub(x).rec)); 
    if isempty(ib); continue; end
    
    st = sub(x).st;
    st_mov = extractEventST(st, beh(ib).on, beh(ib).off, 1); 
    st_imm = extractEventST(st, beh(ib).onRest, beh(ib).offRest, 1); 
    acc = getAcc(beh(ib).vel);
    
    [align{x,1}, time] = getSTA(acc, st_imm, Fs, [win(1), win(2)]);
    align{x,2} = getSTA(acc, st_mov, Fs, [win(1), win(2)]);
    waitbar(x/length(sub),h);
end
fprintf('Done! \n'); close(h);
%     ev_new = shiftST(st, nShuff, 1/nShuff); %shuffle spike times n times
%     mat_shuff = [];
%     for z = 1:nShuff
%         tmp_shuff = getSTA(sig, ev_new{z}, Fs, [time(1), time(end)]);
%         mat_shuff(:,z) = nanmean(tmp_shuff,2); 
%     end
%     mu = nanmean(nanmean(mat_shuff,2)); sigma = nanmean(nanstd(mat_shuff,[],2));
%     mat_z = (mat - nanmean(mu))./nanmean(sigma);
%     align{x} = mat; alignZ{x} = matZ;

align_mov = []; align_imm = [];
for x = 1:length(align)
    if isempty(align{x,2}); align_mov(:,x) = nan(length(time),1); continue; end
    align_mov(:,x) = nanmean(align{x,2},2); 
    if isempty(align{x,1}); align_imm(:,x) = nan(length(time),1); continue; end
    align_imm(:,x) = nanmean(align{x,1},2); 
end

%%
figure; hold on
shadederrbar(time, nanmean(align_imm,2), SEM(align_imm,2), 'k'); 
shadederrbar(time, nanmean(align_mov,2), SEM(align_mov,2), 'g'); 
ylabel('Acceleration'); xlabel('Latency (s)'); xlim([-0.5 0.5]);
title(sprintf('Acc to CIN st (n = %d units)',size(align_mov,2)));

%%
a = max(align_mov([26:76],:));
[~,b] = sort(a);
c = align_mov([26:76],b);

figure;
h = heatmap(c');
h.Title = 'acc to CIN st mov';
h.XLabel = 'Latency to spike (s)'; h.YLabel = 'Unit';
h.GridVisible = 'off'; 
h.Colormap = jet;

%%
figure; plm = floor(sqrt(length(align))); pln = ceil(length(align)/plm);
for x = 1:length(align)
    if isempty(align{x,2}); continue; end
    sp(x) = subplot(plm,pln,x); 
    shadederrbar(time, nanmean(align{x,2},2), SEM(align{x,2},2), 'k'); 
    xlabel('Latency to spike (s)'); xlim([-0.5 0.5]);
    title(sprintf('%s-#%d',sub(x).rec,sub(x).n));
end
