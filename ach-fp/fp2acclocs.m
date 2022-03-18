%% Input Variables
beh = behwt(10:12); sub = cinwt(20:29); %DMS recordings
beh = behwt([1:9,13:28,36:40]); sub = cinwt([1:19,30:98]); %DLS recordings

%% Extract Peaks
for x = 1:length(beh)
    vel = beh(x).vel;
    acc = [vel(1); diff(fliplr(movmean(fliplr(movmean(vel,10)),10)))]; 
    [~,locs] = findpeaks(acc,'MinPeakProminence',1,'MinPeakDistance',1);
    locs = beh(x).time(locs);
	beh(x).acc_locs = locs; clc
end; clc

%% Align acceleration to CIN unit spikes
Fs = 50; time = [-1:1/Fs:1]; %CHANGE: window for STA
align = {length(beh),2};
h = waitbar(0, 'STA: photometry to acceleration peaks');
for x = 1:length(beh)
    ev = beh(x).acc_locs; % Aligning signal to acceleration peaks
    if isempty(ev); continue; end
    for y = 1:length(beh(x).FP)
        sig = beh(x).FP{y}; % Photometry
        [mat,~,matZ] = getSTA(sig, ev, Fs, [time(1), time(end)]);
    %     ev_new = shiftST(st, nShuff, 1/nShuff); %shuffle spike times n times
    %     mat_shuff = [];
    %     for z = 1:nShuff
    %         tmp_shuff = getSTA(sig, ev_new{z}, Fs, [time(1), time(end)]);
    %         mat_shuff(:,z) = nanmean(tmp_shuff,2); 
    %     end
    %     mu = nanmean(nanmean(mat_shuff,2)); sigma = nanmean(nanstd(mat_shuff,[],2));
    %     mat_z = (mat - nanmean(mu))./nanmean(sigma);
        align{x,y} = mat;
    end
    waitbar(x/length(beh),h);
end
close(h); fprintf('Done! \n');

%%
figure; plm = floor(sqrt(length(align))); pln = ceil(length(align)/plm);
for x = 1:length(align)
    if isempty(align{x}); continue; end
    sp(x) = subplot(plm,pln,x); 
    shadederrbar(time, nanmean(align{x},2), SEM(align{x},2), 'k'); 
    xlabel('Latency to Acc Peak (s)'); grid on; xlim([-1 1]);
    ylabel(sprintf('%s (% dF/F)',beh(x).FPnames{1}));
    title(sprintf('%s: %s',beh(x).rec,beh(x).FPnames{y}),'Interpreter','none');
end
