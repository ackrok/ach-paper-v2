a = [1:length(modACh)]; rmv = [3 17 20 22 23 30 36:39 41:43 48 50]; a(rmv) = [];
beh = modACh(a);

%% Extract mu, standard deviation
mu = []; med = []; sig = [];
for x = 1:length(beh)
    fp = beh(x).FP{1}; Fs = 50;
    idx_imm = extractEventST([1:length(fp)]', beh(x).onRest, beh(x).offRest, 1);
    idx_mov = extractEventST([1:length(fp)]', beh(x).on, beh(x).off, 1);
    
    mu(x,1) = nanmean(fp(idx_imm)); mu(x,2) = nanmean(fp(idx_mov));
    sig(x,1) = nanstd(fp(idx_imm)); sig(x,2) = nanstd(fp(idx_mov));
end

tmp = {}; for x = 1:length(beh); tmp{x} = strtok(beh(x).rec,'-'); end
uni = unique(tmp); nAn = length(uni);
sig_an = [];
for x = 1:nAn
    ii = find(strcmp(tmp,uni{x}));
    sig_an(x,:) = nanmean(sig(ii,:),1);
end

%% PLOT
figure; hold on
plot(sig_an', '.-', 'Color', [0 0 0 0.1], 'MarkerSize', 20);
errorbar([0.75 2.25], nanmean(sig_an), SEM(sig_an,1), '.g', 'MarkerSize', 20);
xticks([1 2]); xticklabels({'imm','loc'}); 
ylabel('standard deviation'); ylim([0 8]);
[~,p] = ttest(sig_an(:,1),sig_an(:,2));
title(sprintf('p = %1.4f',p)); axis('square');
title('ACh STD: ttest p = 2.0623e-05 (n = 14 mice)'); 