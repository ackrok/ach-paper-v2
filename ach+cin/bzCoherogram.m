%%
good = [1:length(modAChDA)]; rmv = 42; good(rmv) = [];
beh = modAChDA(good);

rmv = [1 11 15 16 17 20 21 24];
beh = allDMS; beh(rmv) = [];

%%
mat = struct;
h = waitbar(0,'coherogram');
for x = 1:length(beh); aa = [4 3];
  
    mat(x).rec = beh(x).rec; % load recording name
    fp_1 = [beh(x).FP{aa(1)} - nanmean(beh(x).FP{aa(1)})];
    fp_2 = [beh(x).FP{aa(2)} - nanmean(beh(x).FP{aa(2)})]; 
    
    if isfield(beh,'reward')
        idx_rew = extractEventST([1:length(fp_1)]', floor(beh(x).reward)-100, floor(beh(x).reward)+100, 1); % identify sample during reward
    else; idx_rew = []; end
    idx_mov = extractEventST([1:length(fp_1)]', beh(x).on, beh(x).off, 1); % identify sample during locomotion
    idx_mov_nonRew = idx_mov(~ismember(idx_mov, idx_rew)); % exclude reward, include locomotion
    idx_imm = extractEventST([1:length(fp_1)]', beh(x).onRest, beh(x).offRest, 1); % identify sample during locomotion
    idx_imm_nonRew = idx_imm(~ismember(idx_imm, idx_rew)); % exclude reward, include rest
    c = cell(3,1); c{1} = idx_imm_nonRew; c{2} = idx_mov_nonRew; c{3} = idx_rew; % index into cell array for ease of iteration
    %%
    for y = 1:length(c)
        if ~isempty(c{y})
            sig1 = fp_1(c{y}); % extract indexes samples 
            sig2 = fp_2; 
            % sig2 = [fp_2(1); diff(fp_2)]; Fs = 50; % first derivative of DA photometry
            sig2 = sig2(c{y}); % extract indexes samples
            [coher,ph,t,f] = bz_MTCoherogram(sig1,sig2,'frequency',50,'range',[0 10],'window',10,'overlap',5,'step',5,'tapers',[3 5],'pad',0);
            mat(x).coher(:,y) = nanmean(coher,2); % collapse time dimension
            mat(x).phase(:,y) = nanmean(ph,2); % collapse time dimension
        end
    end
    %%
    tmp_coher = []; tmp_phase = []; new_2 = fp_2;
    for s = 1:50 % repeat shuffle N times
        % new_2 = circshift(new_2, Fs);
        new_2 = new_2(randperm(length(new_2)));
        [c,p] = bz_MTCoherogram(fp_1,new_2,'frequency',50,'range',[0 10],'window',10,'overlap',5,'step',5,'tapers',[3 5],'pad',0);
        tmp_coher(:,s) = nanmean(c,2); % collapse time dimension
        tmp_phase(:,s) = nanmean(p,2); % collapse time dimension
    end
    mat(x).coher_shuff = prctile(tmp_coher, [5 50 95], 2); % 5th, 50th, 95th percentiles
    mat(x).phase_shuff = prctile(tmp_phase, [5 50 95], 2); % 5th, 50th, 95th percentiles
    
    waitbar(x/length(beh),h);
end
close(h); fprintf('Coherogram Analysis Done !\n');

% Extract IMM MOV REW
coher_state = cell(3,1); phase_state = cell(3,1);
for x = 1:length(mat)
    for y = 1:3
        if y == 3 && size(mat(x).coher,2) < 3
            coher_state{y}(:,x) = nan(103,1); phase_state{y}(:,x) = nan(103,1);
        else
        coher_state{y}(:,x) = mat(x).coher(:,y);
        phase_state{y}(:,x) = mat(x).phase(:,y);
        end
    end
end
%
rec = {}; for x = 1:length(beh); rec{x} = strtok(beh(x).rec,'-'); end
uni = unique(rec); nAn = length(uni);
coher_an = cell(3,1); phase_an = cell(3,1);
coher_shuff = cell(3,1); phase_shuff = cell(3,1);
for x = 1:nAn
    ii = find(strcmp(rec,uni{x})); % index of matching recordings for this animal
    coh_tmp = [mat(ii).coher_shuff];
    ph_tmp = [mat(ii).phase_shuff];
    for y = 1:3
        coher_an{y}(:,x) = nanmean(coher_state{y}(:,ii),2); % average recordings for each animal, for each behavioral state
        phase_an{y}(:,x) = nanmean(phase_state{y}(:,ii),2); % average recordings for each animal, for each behavioral state
        coher_shuff{y}(:,x) = nanmean(coh_tmp(:,[y:3:size(coh_tmp,2)]),2); % average recordings for each animal, for each shuff percentile
        phase_shuff{y}(:,x) = nanmean(ph_tmp(:,[y:3:size(ph_tmp,2)]),2); % average recordings for each animal, for each shuff percentile
    end
end

%% AVERAGE comparing behavioral states
fig = figure; fig.Position(3) = 1000;
subplot(1,2,1); hold on
clr = {'r','g','b'};
for y = 1:3; shadederrbar(f, nanmean(coher_an{y},2), SEM(coher_an{y},2), clr{y}); end
% for y = 1:3; shadederrbar(f, nanmean(coher_state{y},2), SEM(coher_state{y},2), clr{y}); end
plot([0.5 0.5],[0 1],'--k'); plot([4 4],[0 1],'--k');
shadederrbar(f, nanmean(coher_shuff{2},2), nanmean(coher_shuff{3},2) - nanmean(coher_shuff{2},2), 'k');
legend({'imm','','mov','','rew','','0.5Hz','4Hz'});
ylabel('coherence'); ylim([0 1]); yticks([0:0.2:1]);
title('coherence magnitude: ACh/DA'); xlabel('frequency');

subplot(1,2,2); hold on
for y = 1:3; shadederrbar(f, nanmean(rad2deg(phase_an{y}),2), SEM(rad2deg(phase_an{y}),2), clr{y}); end
% for y = 1:3; shadederrbar(f, nanmean(rad2deg(phase_state{y}),2), SEM(rad2deg(phase_state{y}),2), clr{y}); end
plot([0.5 0.5],[-180 180],'--k'); plot([4 4],[-180 180],'--k');
shadederrbar(f, nanmean(rad2deg(phase_shuff{2}),2), nanmean(rad2deg(phase_shuff{3}),2) - nanmean(rad2deg(phase_shuff{2}),2), 'k');
legend({'imm','','mov','','rew','','0.5Hz','4Hz'});
ylabel('degrees'); ylim([-180 180]); yticks([-180:90:180]);
title('coherence phase');

%% AVERAGE heatmap for ONE behavioral state
y = 1;
fig = figure; fig.Position(3) = 1000;
subplot(1,2,1);
a = coher_an{y};
% a = coher_state{y};
m = max(a); [~,b] = sort(m); c = a(:,b);
h(y) = imagesc(t, f, a, [0 1]);
colorbar; colormap(jet(256));
title('coherence ACh/DA'); axis('square')
xlabel('time (s)'); ylabel('frequency (Hz)');

subplot(1,2,2);
% a = rad2deg(-phase_an{y}); c = a(:,b); % flip axes so red colors are more prominent
a = rad2deg(-phase_state{y}); c = a(:,b);
h(y) = imagesc(t, f, a, [-180 180]);
colorbar; colormap(jet(256));
title('phase IMM'); axis('square')
xlabel('time (s)'); ylabel('frequency (Hz)');
movegui(gcf,'center');

%% TESTING
x = 18;

fp_mat = [beh(x).FP{1}, beh(x).FP{2}]; % extract photometry signal from structure
fp_mat = fp_mat - nanmean(fp_mat);
Fs = beh(x).Fs;

sig = fp_mat; r = [1 1800];
% r = [1210 1590]; % x = 28, REW + IMM + MOV
% r = [1 600];  % x = 7,  IMMOBILITY
% r = [1200 1600];  % x = 7, IMMOBILITY + two movement bouts
% r = [1400 1800]; % DMS x = 5 (GZ001 210416) IMM + MOV
r = [300 700]; % DMS x = 18 (GZ002 210222) IMM + MOV
% r = [250 650]; % openfield x = 7 GZ002-210426

sig = fp_mat([find(beh(x).time == r(1)):find(beh(x).time == r(2))],:);

%
[coher,ph,t,f] = bz_MTCoherogram(sig(:,1),sig(:,2),'frequency',50,'range',[0 10],'window',10,'overlap',5,'step',5,'tapers',[3 6],'pad',0);

%
fig = figure; fig.Position([3 4]) = [1200 1200];

subplot(2,2,1);
y = 1; a = coher;
sm = 0.75;
a = imgaussfilt(a, sm); % filter image with Gaussian smoothing kernel of X std
h(y) = imagesc(t+r(1), f, a, [0 1]);
colorbar; colormap(jet(256));
title(sprintf('%s coherence ACh/DA',beh(x).rec)); axis('square')
xlabel('time (s)'); ylabel('frequency (Hz)');

subplot(2,2,2);
y = 2; a = rad2deg(-ph); % flip axes so red colors are more prominent
a = imgaussfilt(a, sm); % filter image with Gaussian smoothing kernel of X std\
h(y) = imagesc(t+r(1), f, a, [-180 180]);
colorbar; colormap(jet(256));
title(sprintf('%s phase IMM',beh(x).rec)); axis('square')
xlabel('time (s)'); ylabel('frequency (Hz)');

subplot(2,2,3); hold on
r = [find(beh(x).time == r(1)) find(beh(x).time == r(2))];
plot(beh(x).time(r(1):r(2)), getAcc(beh(x).vel(r(1):r(2))), 'k');
on = beh(x).on(beh(x).on > r(1) & beh(x).on < r(2));
off = beh(x).off(beh(x).off > r(1) & beh(x).off < r(2));
stem(on./Fs, ones(length(on),1), 'g');
stem(off./Fs, ones(length(off),1), 'r');
if ~isempty(beh(x).reward); rew = beh(x).reward(beh(x).reward > r(1) & beh(x).reward < r(2)); 
    stem(rew./Fs, ones(length(rew),1), 'b'); end
ylabel('Acceleration'); ylim([-1 1]); axis('square');

subplot(2,2,4); hold on
plot(beh(x).time(r(1):r(2)), sig(:,1), 'g');
plot(beh(x).time(r(1):r(2)), sig(:,2), 'm');
ylabel('Fluorescence'); axis('square');

movegui(gcf, 'center');

%%
fig = figure; fig.Position(3) = 1200;
subplot(1,2,1);
y = 1; h(y) = pcolor(t, f, coher);
shading interp
set(h(y),'EdgeColor','none','linestyle','none');
colorbar; colormap(jet(256));
title(sprintf('%s coherence magnitude ACh/DA',beh(x).rec)); axis('square')
xlabel('time (s)'); ylabel('frequency (Hz)');

subplot(1,2,2);
y = 2; h(y) = heatmap(rad2deg(ph));
y = 1; h(y) = pcolor(t, f, rad2deg(ph));
shading interp
set(h(y),'EdgeColor','none','linestyle','none');
colorbar; colormap(jet(256));
title(sprintf('%s coherence magnitude ACh/DA',beh(x).rec)); axis('square')
xlabel('time (s)'); ylabel('frequency (Hz)');

%%
figure;
yyaxis left; plot(f, nanmean(coher,2)); ylabel('coherence'); ylim([0 1]);
yyaxis right; plot(f, rad2deg(nanmean(ph,2))); ylabel('degrees'); ylim([-180 180]);
xlabel('frequency (Hz)');
%%
