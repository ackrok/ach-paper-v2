figure; plm = floor(sqrt(length(beh))); pln = ceil(length(beh)/plm);
for x = 1:length(beh)
reward = beh(x).reward; lick = beh(x).lick; Fs = beh(x).Fs;
bin = 0.1; window = [-2 2];
peth = getClusterPETH(lick/Fs, reward/Fs, bin, window);
sp(x) = subplot(plm,pln,x);
shadederrbar(peth.time, nanmean(peth.cts{1},2), SEM(peth.cts{1},2), 'm');
title(sprintf('%s',beh(x).rec));
end
subplot(plm,pln,1);
xlabel('Latency to reward (s)'); ylabel('Lick');