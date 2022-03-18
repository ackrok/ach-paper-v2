%% (1) Load data
load('fullcin+beh_ACh')
sub = cinACh;  % Extract CINs from recordings with ACh
beh = behACh; % Extract ACh recordings

%% (2) PETH
% bin = 0.05; window = [-0.025 0.025]; %CHANGE: window for PETH
bin = 0.02; window = [-0.01 0.01]; %CHANGE: window for PETH

uni = unique({sub.rec}); % How many unique recordings are there?
mat = struct; % Initialize structure
h = waitbar(0, 'PETH: CIN spikes to FP peak times');
for u = 1:length(uni)
    ii = find(strcmp({sub.rec},uni{u})); % Index of units from this recording
    ib = find(strcmp({beh.rec},uni{u})); % Index of matching photometry data
    if length(ii) < 2; continue; end % Skip recordings where less than 2 units
    sub_uni = sub(ii); % New structure with only units from this recording
    h2 = waitbar(0, 'PETH: this recording');
    fprintf('%s ... ',uni{u});
    for x = 1:length(ii)
        %% Extract spike times
        jj = [1:length(ii)]; jj(x) = []; % "other" 
        st = [sub_uni(x).st]; % Extract spike times of this unit
        st_other = {sub_uni(jj).st}; % Extract spike times of all other units from this recording
        %% PETH
        peth = getClusterPETH(st_other, st, bin, window); % PETH: spike times aligned to spike times
        %%
        prob = []; 
        t_0 = 1; cts0 = []; cts = {};
        for y = 1:length(st_other)
            cts{y} = peth.cts{y}; 
            prob(:,y) = (sum(peth.cts{y},2))./length(st_other{y}); % Probability of firing = cts/st
            cts0(y,:) = peth.cts{y}(t_0,:);
        end 
        %% STA
        fp = beh(ib).FP{1}; Fs = 50;
        fp = fp - nanmean(fp);
        % fp = [beh(1).vel(1); diff(movmean(beh(1).vel,10))]; Fs = 50;
        [sta_fp, sta_time] = getSTA(fp, st, Fs, [-1, 2]);
        
        %% Load into output structure
        next = 1 + length(mat);
        mat(next).rec = sub_uni(x).rec; 
        mat(next).n = sub_uni(x).n; mat(next).m = [sub_uni(jj).n];
        mat(next).st = st;
        mat(next).cts0 = cts0; 
        mat(next).cts = cts; mat(next).prob = prob; 
        mat(next).sta = sta_fp;
        
        %% REST vs MVMT
        st_imm = extractEventST(st, beh(ib).onRest, beh(ib).offRest, 1);
        idx_imm = find(ismember(st, st_imm)); % Idx of spike times during immobility
        st_mov = extractEventST(st, beh(ib).on, beh(ib).off, 1);
        idx_mov = find(ismember(st, st_mov)); % Idx of spike times during locomotion
        mat(next).idx_imm = idx_imm; mat(next).idx_mov = idx_mov;
        
        %%
        waitbar(x/length(ii),h2);
    end; close(h2);
    fprintf('done. \n');
    waitbar(u/length(uni),h);
end; close(h);
if isempty(mat(1).n); mat(1) = []; end
fprintf('ANALYSIS COMPLETE: aligning spike times to spike times \n');
time = peth.time; 

%% SHUFF
nShuff = 10; % change number of shuffles
for x = 1:length(mat)
    st = mat(x).st;
    st_imm = st(mat(x).idx_imm); if isempty(st_imm); continue; end
    shuffSt = shuffleST(st_imm, nShuff); % shuff spike times
    ib = find(strcmp({beh.rec},mat(x).rec));
    fp = beh(ib).FP{1}; Fs = 50;
    fp = fp - nanmean(fp);
    prc5 = []; prc50 = []; prc95 = [];
    for y = 1:nShuff
        tmp2 = getSTA(fp, shuffSt{y}, Fs, [-1, 2]); % align photometry to shuffled spike times
        prc = prctile(tmp2, [5 50 95], 2);
        prc5(:,y) = prc(:,1);
        prc50(:,y) = prc(:,2);
        prc95(:,y) = prc(:,3);
    end
    tmp = [nanmean(prc5,2), nanmean(prc50,2), nanmean(prc95,2)];
    mat(x).prc_shuffRest = tmp./nShuff; 
end
fprintf('SHUFFLE DONE. \n');

%% (3) 1/N vs N/N IMMOBILITY
sta_max = cell(length(mat),2);
sta_max_imm = nan(size(mat(1).sta,1),length(mat)); sta_zero_imm = sta_max_imm; 
sta_shuff = sta_max_imm; sta_shuff95 = sta_max_imm;
propN = [];
for x = 1:length(mat)
    if isempty(mat(x).idx_imm); continue; end
    cts = mat(x).cts0; % extract counts
    cts(cts > 1) = 1; 
    cts = sum(cts, 1); % sum across units
    iN = find(cts >= length(mat(x).m)); iN = iN'; % index of max coherence among units
    i0 = find(cts == 0); i0 = i0'; % index of no coherence
    iN_imm = mat(x).idx_imm((ismember(mat(x).idx_imm, iN))); % index of max coherence among units that occur during immobility
    i0_imm = mat(x).idx_imm((ismember(mat(x).idx_imm, i0))); % index of zero coherence among units that occur during immobility
    pull_sta = mat(x).sta;
    pull_sta = pull_sta - nanmean(pull_sta([1:10],:));
    sta_max{x,1} = pull_sta(:,i0_imm);
    sta_max{x,2} = pull_sta(:,iN_imm);
    sta_zero_imm(:,x) = nanmean(sta_max{x,1},2);
    sta_max_imm(:,x) = nanmean(sta_max{x,2},2);
    sta_shuff(:,x) = mat(x).prc_shuffRest(:,2);
    sta_shuff95(:,x) = mat(x).prc_shuffRest(:,1) - mat(x).prc_shuffRest(:,2);
    propN(x) = length(iN_imm)/length(cts);
end

%% (3b) Plot STA IMM for each recording - IMMOBILITY
figure; plm = floor(sqrt(length(mat))); pln = ceil(length(mat)/plm);
for x = 1:length(mat)
    if isempty(mat(x).idx_imm); continue; end
    sp(x) = subplot(plm, pln, x); hold on
    shadederrbar(sta_time, mat(x).prc_shuffRest(:,2), mat(x).prc_shuffRest(:,1) - mat(x).prc_shuffRest(:,2),'k'); 
    shadederrbar(sta_time, nanmean(sta_max{x,1},2), SEM(sta_max{x,1},2),'m');
    shadederrbar(sta_time, nanmean(sta_max{x,2},2), SEM(sta_max{x,2},2),'g');
    xlim([-1 1]);
    title(sprintf('%s-%d',mat(x).rec,mat(x).n),'Interpreter','none');
end
linkaxes(sp,'y');

%% (3c) N/N average across all recordings - IMMOBILITY
% figure; hold on
fig = figure; fig.Position(3) = 1000;
subplot(1,2,1); hold on
plot([0 0],[-0.5 1.0],'--k'); 
shadederrbar(sta_time, nanmean(sta_shuff,2), nanmean(sta_shuff95,2), 'k'); hold on
shadederrbar(sta_time, nanmean(sta_zero_imm,2), SEM(sta_zero_imm,2), 'm'); hold on
shadederrbar(sta_time, nanmean(sta_max_imm,2), SEM(sta_max_imm,2), 'g');
ylabel('ACh3.0 (delta %F/F)'); ylim([-0.5 1.5]); 
xlabel('Latency to CIN spike (s)'); xlim([-1 1]);
title(sprintf('ACh to 1/N vs N/N CIN spikes (n = %d units)',size(sta_max_imm,2)));
axis('square');

%% (3d) how many >95% CI - IMMOBILITY
above95_max = []; above95_zero = []; below5_max = []; below5_zero = [];
for x = 1:length(mat)
    if isempty(mat(x).idx_imm); continue; end
    above95_zero = [above95_zero, sta_zero_imm(:,x) > mat(x).prc_shuffRest(:,3)];
    above95_max  = [above95_max,  sta_max_imm(:,x)  > mat(x).prc_shuffRest(:,3)];
    below5_zero  = [below5_zero,  sta_zero_imm(:,x) < mat(x).prc_shuffRest(:,1)];
    below5_max   = [below5_max,   sta_max_imm(:,x)  < mat(x).prc_shuffRest(:,1)];
end
ds = 2;
above95_max = above95_max(1:ds:end,:); above95_zero = above95_zero(1:ds:end,:); 
below5_max = below5_max(1:ds:end,:); below5_zero = below5_zero(1:ds:end,:);

% figure; hold on
subplot(1,2,2); hold on
bar(sta_time(1:ds:end), 100*sum(above95_zero,2)/size(above95_zero,2),'FaceColor','m','FaceAlpha',0.5);
bar(sta_time(1:ds:end), -100*sum(below5_zero,2)/size(below5_zero,2),'FaceColor','m','FaceAlpha',0.5);
bar(sta_time(1:ds:end), 100*sum(above95_max,2)/size(above95_max,2),'FaceColor','g','FaceAlpha',0.5);
bar(sta_time(1:ds:end), -100*sum(below5_max,2)/size(below5_max,2),'FaceColor','g','FaceAlpha',0.5);
xlabel('Latency to CIN spike (s)'); xlim([-1 1]);
ylabel('Prop > 95% CI (%)'); ylim([-100 100]);
title(sprintf('Proportion > 95p CI (n = %d units) - IMM',size(below5_max,2)))
axis('square');

%% (3e) EXAMPLE
x = 4;
figure; hold on
% subplot(1,3,1); hold on
clr = {'m','r','c','b','g'};
cts = mat(x).cts0; % extract counts
cts(cts > 1) = 1; 
cts = sum(cts, 1); % sum across units
for y = 1:length(find(strcmp({mat.rec},mat(x).rec)))
    iTmp = find(cts == y-1); iTmp = iTmp'; % index of coherence among units
    iTmp_imm = mat(x).idx_imm((ismember(mat(x).idx_imm, iTmp))); % index of coherence among units that occur during immobility
    pull_sta = mat(x).sta;
    pull_sta = pull_sta - nanmean(pull_sta([1:10],:));
    
    shadederrbar(sta_time, mat(x).prc_shuffRest(:,2), mat(x).prc_shuffRest(:,1) - mat(x).prc_shuffRest(:,2), 'k');
    shadederrbar(sta_time, nanmean(pull_sta(:,iTmp_imm),2), SEM(pull_sta(:,iTmp_imm),2), clr{y});
end
xlabel('Latency to CIN spike (s)'); xlim([-1 1]); 
ylabel('Photometry (% dF/F)');
title(sprintf('%s-#%d',mat(x).rec,mat(x).n),'Interpreter','none');
axis('square');
movegui(fig, 'center');

%% (3f) DA latency to pCIN spike
max_lat = nan(length(mat),1); max_val = max_lat;
for x = 1:length(mat)
    if isempty(mat(x).idx_imm); continue; end
    [mm, ii] = min(nanmean(sta_max{x,2},2)); % min of N/N coherence
    max_lat(x) = sta_time(ii);
    max_val(x) = mm;
end
figure; violinplot(max_lat);

%% (4) 1/N vs N/N -- LOCOMOTION
% sta_max = cell(length(mat),2);
% sta_max_mov = []; sta_zero_mov = []; 
% for x = 1:length(mat)
%     cts = mat(x).cts0; % extract counts
%     cts(cts > 1) = 1; 
%     cts = sum(cts, 1); % sum across units
%     iN = find(cts >= length(mat(x).m)); iN = iN'; % index of max coherence among units
%     i0 = find(cts == 0); i0 = i0'; % index of no coherence
%     iN_mov = mat(x).idx_mov((ismember(mat(x).idx_mov, iN))); % index of max coherence among units that occur during immobility
%     i0_mov = mat(x).idx_mov((ismember(mat(x).idx_mov, i0))); % index of zero coherence among units that occur during immobility
%     pull_sta = mat(x).sta;
%     pull_sta = pull_sta - nanmean(pull_sta([1:10],:));
%     sta_max{x,1} = pull_sta(:,i0_mov);
%     sta_max{x,2} = pull_sta(:,iN_mov);
%     sta_zero_mov(:,x) = nanmean(sta_max{x,1},2);
%     sta_max_mov(:,x) = nanmean(sta_max{x,2},2);
% end
% 
% figure; hold on
% plot([0 0],[-0.5 1.0],'--k'); 
% shadederrbar(sta_time, nanmean(sta_shuff,2), nanmean(sta_shuff95,2), 'k'); hold on
% shadederrbar(sta_time, nanmean(sta_zero_mov,2), SEM(sta_zero_mov,2), 'm'); hold on
% shadederrbar(sta_time, nanmean(sta_max_mov,2), SEM(sta_max_mov,2), 'g');
% ylabel('ACh3.0 (delta %F/F)'); ylim([-0.5 1.0]); 
% xlabel('Latency to CIN spike (s)');
% title(sprintf('ACh to 1/N vs N/N CIN spikes (n = %d units) - MOV',size(sta_max_mov,2)));
% 
% above95_max = []; above95_zero = []; below5_max = []; below5_zero = [];
% for x = 1:length(mat)
%     above95_zero = [above95_zero, sta_zero_mov(:,x) > mat(x).prc_shuff(:,3)];
%     above95_max  = [above95_max,  sta_max_mov(:,x)  > mat(x).prc_shuff(:,3)];
%     below5_zero  = [below5_zero,  sta_zero_mov(:,x) < mat(x).prc_shuff(:,1)];
%     below5_max   = [below5_max,   sta_max_mov(:,x)  < mat(x).prc_shuff(:,1)];
% end
% figure; hold on
% bar(sta_time, 100*sum(above95_zero,2)/size(above95_zero,2),'FaceColor','m','FaceAlpha',0.5);
% bar(sta_time, -100*sum(below5_zero,2)/size(below5_zero,2),'FaceColor','m','FaceAlpha',0.5);
% bar(sta_time, 100*sum(above95_max,2)/size(above95_max,2),'FaceColor','g','FaceAlpha',0.5);
% bar(sta_time, -100*sum(below5_max,2)/size(below5_max,2),'FaceColor','g','FaceAlpha',0.5);
% xlabel('Latency to CIN spike (s)'); xlim([-1 1]);
% ylabel('Prop > 95% CI (%)'); ylim([-100 100]);
% title(sprintf('Proportion > 95p CI (n = %d units) - MOV',size(below5_max,2)))