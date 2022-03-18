clr = [0.05 0.75 0.45; 0.5 0.2 0.55; 0.85 0.35 0.1; 0.3 0.75 0.9];

%% PLOT - by animal
z = 1; % 1: pause, 2: peak
switch z; case 1; lbl = 'Pause'; case 2; lbl = 'Peak'; end 

fig = figure; fig.Position([3 4]) = [600 600];
for x = 1:4
    subplot(2,2,x);
    for y = [1 2 4]
        a = s(y).da2ach{x,z}; % Extract aligned photometry signal
        % a = s(y).ach2ach{x,z};
        if isempty(a); continue; end
        a = a - nanmean(a([1:find(sta_time == -1)],:)); % Subtract baseline from -2:-1
        shadederrbar(sta_time, nanmean(a,2), SEM(a,2), clr(y,:));
    end
    xlabel('Latency (s)'); xlim([-1 1]);
    ylabel('rDA1m (%dF/F)'); title(sprintf('DA to ACh %s',lbl)); axis('square');
    % ylabel('ACh3.0 (%dF/F)'); title(sprintf('ACh to ACh %s',lbl)); axis('square');
end
movegui('center');

%% AVERAGE
nAn = 4; 
da2ach_avg = cell(4,2); 
da2ach_max = cell(4,2);
time_sub = sta_time(find(sta_time == -1):find(sta_time == 1));
for z = 1:2
    for y = 1:4
        for x = 1:nAn
            a = s(y).da2ach{x,z}; % Extract aligned photometry signal
            if isempty(a)
                da2ach_avg{y,z}(:,x) = nan(101,1); 
                da2ach_max{y,z}(x) = nan;
                continue; end
            a = a - nanmean(a([1:find(sta_time == -1)],:)); % Subtract baseline from -2:-1
            a = a(find(sta_time == -1):find(sta_time == 1),:);
            da2ach_avg{y,z}(:,x) = nanmean(a,2); % Store average across all trials
            da2ach_max{y,z}(x) = max(abs(da2ach_avg{y,z}([find(time_sub == -0.5):find(time_sub == 0)],x)));
        end
    end
end
time_sub = sta_time(find(sta_time == -1):find(sta_time == 1));

%% PLOT AVERAGE
fig = figure; fig.Position(3) = 1000;
for z = 1:2
    switch z; case 1; lbl = 'Pause'; case 2; lbl = 'Peak'; end 
    subplot(1,2,z); hold on
    for y = [1 2 4]
        shadederrbar(time_sub, nanmean(da2ach_avg{y,z},2), SEM(da2ach_avg{y,z},2), clr(y,:));
    end
    xlabel('Latency (s)'); xlim([-1 1]);
    ylabel('rDA1m (%dF/F)');
    title(sprintf('DA to ACh %s',lbl));
end

%% PLOT MAXIMUM
fig = figure; fig.Position(3) = 1000;
for z = 1:2
    switch z; case 1; lbl = 'Pause'; case 2; lbl = 'Peak'; end 
    subplot(1,2,z); hold on
    a = [da2ach_max{1,z}',da2ach_max{2,z}',da2ach_max{4,z}'];
    a = a./nanmean(a(:,1));
    
    plot([1 3],[1 1],'--k');
    errorbar(abs(nanmean(a)'),abs(SEM(a,1)'),'.k','MarkerSize',20);
    plot([1.15; 1.85; 2.85].*ones(3,nAn),abs(a'),'.k','MarkerSize',20);
    plot([1.15; 1.85].*ones(2,nAn),abs(a(:,[1 2])'),':m','MarkerSize',20);
    plot([1.15; 2.85].*ones(2,nAn),abs(a(:,[1 3])'),':b','MarkerSize',20);
    xlim([0.5 3.5]); xticks([1:3]); xticklabels({'aCSF','D1R/D2R','DHBE'});
        
    p = []; [~,p(2)] = ttest(a(:,1),a(:,2)); [~,p(3)] = ttest(a(:,1),a(:,3));
    ylabel(sprintf('rDA1m to ACh %s',lbl)); ylim([0 4]); 
    title(sprintf('DA to ACh %s (p = %1.4f, %1.4f)',lbl,p(2),p(3)));
    axis('square');
end