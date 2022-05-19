mat = struct;
corr_cell = cell(3,4);
for a = 1:3; for b = 1:4; corr_cell{a,b} = nan(501,length(beh)); end; end

h = waitbar(0, 'cross corr');
for x = 1:length(beh); y = [2 1]; %CHANGE - which FP signal to run CCG or ACG on
    
    %% extract signals
    fp_mat = [];
    fp_mat(:,1) = beh(x).FP{y(1)}; Fs = beh(x).Fs; % extract photometry signal from structure
    fp_mat(:,1) = fp_mat(:,1) - nanmean(fp_mat(:,1)); % subtract baseline (mean of entire photometry signal) from fp
    fp_mat(:,2) = beh(x).FP{y(2)}; % dopamine
    fp_mat(:,2) = fp_mat(:,2) - nanmean(fp_mat(:,2)); % subtract baseline (mean of entire photometry signal) from fp

    %%
    mat(x).rec = beh(x).rec;
      
    [corr_tmp, lags] = xcorr(fp_mat(:,1), fp_mat(:,2), 5*Fs, 'coeff'); % cross-correlation
        
    fp_sub_new = fp_mat(:,2);
    tmp_shuff = []; 
    for s = 1:50
        fp_sub_new = circshift(fp_sub_new, Fs);
        % tmp_shuff(:,s) = xcorr(fp_sub(randperm(size(fp_sub,1)),1), fp_sub(randperm(size(fp_sub,2)),2), 10*Fs, 'coeff');
        % tmp_shuff(:,s) = xcorr(fp_sub(:,1), fp_sub(randperm(size(fp_sub,2)),2), 10*Fs, 'coeff');
        tmp_shuff(:,s) = xcorr(fp_mat(:,1), fp_sub_new, 5*Fs, 'coeff');
    end
        
    mat(x).corr = corr_tmp;
    mat(x).shuff = prctile(tmp_shuff, [5 50 95], 2);
    corr_cell{1}(:,x) = corr_tmp;       % cross-correlation
    corr_cell{2}(:,x) = prctile(tmp_shuff, 5, 2); % shuffle 5th percentile
    corr_cell{3}(:,x) = prctile(tmp_shuff, 50, 2); % shuffle 50th percentile
    corr_cell{4}(:,x) = prctile(tmp_shuff, 95, 2); % shuffle 95th percentile
    
%%
    waitbar(x/length(beh),h);
end
close(h);

%N = X mice
corr_an = cell(3,4); min_val = []; min_lag = [];
tmp = {}; for x = 1:length(beh); tmp{x} = strtok(beh(x).rec,'-'); end
uni = unique(tmp); nAn = length(uni);

for x = 1:nAn
    idx = strcmp(tmp,uni{x});
    for b = 1:4
        corr_adj = corr_cell{b};
        if b == 1; corr_adj = corr_adj - nanmean(corr_adj([1:find(lags == -2)],:)); end
        corr_an{b}(:,x) = nanmean(corr_adj(:,idx),2);
    end
end

[min_val, ii] = min(corr_an{1});
min_lag(:,z) = lags(ii)./50;

%%
mat = struct;
h = waitbar(0,'coherogram');
for x = 1:length(beh); aa = [2 1];
  
    mat(x).rec = beh(x).rec; % load recording name
    fp_1 = [beh(x).FP{aa(1)} - nanmean(beh(x).FP{aa(1)})];
    fp_2 = [beh(x).FP{aa(2)} - nanmean(beh(x).FP{aa(2)})]; 
    
    %%
    [coher,ph,t,f] = bz_MTCoherogram(fp_1,fp_2,'frequency',50,'range',[0 10],'window',10,'overlap',5,'step',5,'tapers',[3 5],'pad',0);
    mat(x).coher = nanmean(coher,2); % collapse time dimension
    mat(x).phase = nanmean(ph,2); % collapse time dimension

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
coher_state = cell(1); phase_state = cell(1);
for x = 1:length(mat); y = 1;
    coher_state{y}(:,x) = mat(x).coher;
    phase_state{y}(:,x) = mat(x).phase;
end
%
rec = {}; for x = 1:length(beh); rec{x} = strtok(beh(x).rec,'-'); end
uni = unique(rec); nAn = length(uni);
coher_an = cell(1); phase_an = cell(1);
coher_shuff = cell(1); phase_shuff = cell(1);
for x = 1:nAn
    ii = find(strcmp(rec,uni{x})); % index of matching recordings for this animal
    coh_tmp = [mat(ii).coher_shuff];
    ph_tmp = [mat(ii).phase_shuff];
    y = 1;
    coher_an{y}(:,x) = nanmean(coher_state{y}(:,ii),2); % average recordings for each animal, for each behavioral state
    phase_an{y}(:,x) = nanmean(phase_state{y}(:,ii),2); % average recordings for each animal, for each behavioral state
    coher_shuff{y}(:,x) = nanmean(coh_tmp(:,[y:3:size(coh_tmp,2)]),2); % average recordings for each animal, for each shuff percentile
    phase_shuff{y}(:,x) = nanmean(ph_tmp(:,[y:3:size(ph_tmp,2)]),2); % average recordings for each animal, for each shuff percentile
end

%%
fig = figure; fig.Position(3) = 1375;
subplot(1,3,1); hold on
plot(lags/50, corr_an{1}, 'Color', [0 0 0 0.2]); 
shadederrbar(lags/50, nanmean(corr_an{1},2), SEM(corr_an{1},2), 'b'); 
% a = nanmean(corr_an{1,3},2);
% a = a - nanmean(a(find(lags/Fs == -5):find(lags/Fs == -1),:));
% shadederrbar(lags/Fs, a, nanmean(corr_an{1,4}-corr_an{1,3},2), 'k'); hold on
xlim([-1 1]); ylim([-0.5 0.25]); yticks([-0.5:0.25:0.5]);
title('coefficient'); axis('square');

subplot(1,3,2); hold on
plot(f, coher_an{1}, 'Color', [0 0 0 0.2]); 
shadederrbar(f, nanmean(coher_an{1},2), SEM(coher_an{1},2), 'b'); 
shadederrbar(f, mat(1).coher_shuff(:,2), mat(1).coher_shuff(:,3) - mat(1).coher_shuff(:,2), 'k');
% shadederrbar(f, nanmean(coher_shuff{1},2), SEM(coher_shuff{1},2), 'k'); 
plot([0.5 0.5],[0 1],'--k'); plot([4 4],[0 1],'--k');
ylim([0 1]);
title('coherence'); axis('square');

subplot(1,3,3); hold on
plot(f, rad2deg(phase_an{1}), 'Color', [0 0 0 0.2]); 
shadederrbar(f, nanmean(rad2deg(phase_an{1}),2), SEM(rad2deg(phase_an{1}),2), 'b'); 
shadederrbar(f, rad2deg(mat(1).phase_shuff(:,2)), rad2deg(mat(1).phase_shuff(:,3) - mat(1).phase_shuff(:,2)), 'k');
% shadederrbar(f, nanmean(rad2deg(phase_shuff{1}),2), SEM(rad2deg(phase_shuff{1}),2), 'k'); 
plot([0.5 0.5],[-180 180],'--k'); plot([4 4],[-180 180],'--k');
ylim([-180 180]); yticks([-180:90:180]);
title('phase'); axis('square');

movegui(gcf,'center');