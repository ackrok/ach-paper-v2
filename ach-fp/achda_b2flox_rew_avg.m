beh = b2flox([19 22 24:26 28:30]); 
[align_full] = plot_fp2event(beh,[-6 2],0);

%% n = X mice, concatenate all trials for each animal
an = {}; for x = 1:length(beh); an{x} = strtok(beh(x).rec,'-'); end
uni = unique(an);

align_adj = {}; % adjust for baseline
for x = 1:length(beh); for y = 1:2; align_adj{x,y} = align_full{x,y} - nanmean(beh(x).FP{y}); end; end

align_an = {}; % Cell array with concatenated trials for each animal
for x = 1:length(uni) % Iterate over each unique animal
    for y = 1:2 % Repeat for ACh, DA
        align_an{x,y} = [];
        matchAN = find(strcmp(an,uni{x})); % Find all recordings that match animal ID
        for z = 1:length(matchAN)
            align_an{x,y} = [align_an{x,y}, align_adj{matchAN(z),y}]; % Concatenate trials from recordings with same animal ID
        end
    end
end

% n = X mice, average all trials for each animals
align_an_avg = cell(1,2); 
align_an_avg_norm = cell(1,2);
for x = 1:size(align_an,1)
    for y = 1:size(align_an,2)
        tmp = align_an{x,y};
        align_an_avg{1,y}(:,x) = nanmean(tmp,2); % Average over all trials for each animal
        tmp = nanmean(tmp,2);
        bb = nanmean(tmp([find(t==-6):find(t==-1)],:),1);
        
        if y == 1 
            mm = min(tmp([find(t==0):find(t==1)],:));
            norm = (tmp - mm)./(bb - mm);
        elseif y == 2
            mm = max(tmp([find(t==0):find(t==1)],:));
            norm = (tmp - bb)./(mm - bb); end
        align_an_avg_norm{1,y}(:,x) = norm; % normalize(align_an_avg{1,y}(:,x),'range'); % Normalize so range is [0 1]
    end
end

%% PLOT: %dF/F 
fig = figure; fig.Position(3) = 1375;
subplot(1,3,1); hold on
plot([0 0],[-5 11],'--k'); hold on
shadederrbar(t, nanmean(align_an_avg{2},2), SEM(align_an_avg{2},2), 'm'); hold on
shadederrbar(t, nanmean(align_an_avg{1},2), SEM(align_an_avg{1},2), 'g');
ylabel('FP (% dF/F)'); xlabel('Latency to reward delivery (s)');
xlim([-1 2]);
axis('square')
title(sprintf('reward: ACh3.0 + rDA1m (n = %d mice)',size(align_an,1)));
 