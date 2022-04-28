% Assuming you have all the pre-infusion values for each drug condition, 
% I would try to normalize the infusion values for each drug to the mode of 
% their respective baseline (pre-infusion) to get rid of differences in 
% expression levels/light levels between animals. Only then would you be 
% able to compare acsf vs. other antagonists.
% 
% load('C:\Users\Anya\Desktop\FP_LOCAL\AK189-197_cannula+rew_v2.mat')
% load('C:\Users\Anya\Desktop\FP_LOCAL\cannula\scop_processed_mod.mat')

a = [34 36 38 40; 26 42 30 44; 18 20 22 24; 46 48 NaN 50]';
tmp = []; tmp_base = [];

for x = 1:4; for y = 1:4
    %%
    if isnan(a(x,y)); tmp(x,y) = nan; tmp_base(x,y) = nan; continue; end 
    demod = cannula(a(x,y)).nbFP{1}; % Extract demod (demodulated, non-baselined) photometry signal for cannula recording
    % demod = cannula(a(x,y)).FP{1}; 
    idxImm = extractEventST([1:length(demod)]',cannula(a(x,y)).onRest,cannula(a(x,y)).offRest,1); % Sample index during immobility
    idxImm(idxImm < s(y).win(x,1)*60*50 | idxImm > s(y).win(x,2)*60*50) = []; % Remove sample indices that are outside of infusion window
    demodImm = demod(idxImm); % Restrict to immobility + infusion window
    base = cannula(a(x,y)-1).nbFP{1}; % Extract demod vector for baseline recording
    % base = cannula(a(x,y)-1).FP{1};
    baseImm = base(extractEventST([1:length(base)]',cannula(a(x,y)-1).onRest,cannula(a(x,y)-1).offRest,1)); % Restrict to immobility
    
    %%
    tmp(x,y) = mode(demodImm); %([1:5*50])); % Mean of first 5s of immobility
    tmp_base(x,y) = mode(baseImm); %([1:5*50])); 
    end; end
%%
scop = [];
for x = 1:6
    demod = scopmod(x).nbFP{1};
    % demod = scopmod(x).FP{1};
    idxImm = extractEventST([1:length(demod)]',scopmod(x).onRest,scopmod(x).offRest,1); % Sample index during immobility
    scop(x) = mode(demod(idxImm));
end
scop_per = scop(4:6)./scop(1:3);

%% FRACTION of BASELINE
figure; hold on
b = [[tmp./tmp_base], [scop_per(:);nan]];
plot([1:5]'.*ones(5,4), b', '.k', 'MarkerSize', 20);
errorbar(0.25+[1:5],nanmean(b,1),SEM(b,1),'.b', 'MarkerSize', 20);
xlim([0.5 5.5]); xticks([1:5]); 
xticklabels({'aCSF','d1d2','glu','nAChR','mAChRantag'});
ylabel('% of baseline'); ylim([0 1]); yticks([0:0.25:1]);
title('mode(infusion)/mode(baseline) IMM');

%% FRACTION of BASELINE -- to aCSF
figure; hold on
b = tmp./tmp_base; b = b./nanmean(b(:,1));
plot([1 4],[1 1],'--k');
plot([1;2;3;4].*ones(4,4), b', '.k', 'MarkerSize', 20);
errorbar([1 2 3 4],nanmean(b,1),SEM(b,1),'b');
xlim([0.5 4.5]); xticks([1:4]); xticklabels({'aCSF','d1d2','glu','nAChR'});
ylabel('% of baseline (norm.)'); ylim([0 2]); yticks([0:0.25:2]);
title('mode(infusion)/mode(baseline) IMM');

%% FRACTION of BASELINE -- 1 set to aCSF and 0 set to mAChR antagonist
figure; hold on
b = tmp./tmp_base;
b = (b - nanmean(scop))./(nanmean(tmp(:,1)./tmp_base(:,1)) - nanmean(scop));
plot([1 4],[1 1],'--k');
plot([1;2;3;4].*ones(4,4), b', '.k', 'MarkerSize', 20);
errorbar([1 2 3 4],nanmean(b,1),SEM(b,1),'b');
xlim([0.5 4.5]); xticks([1:4]); xticklabels({'aCSF','d1d2','glu','nAChR'});
ylabel('% of baseline (norm. 1-aCSF 0-mAChR)'); ylim([0 1.5]); yticks([0:0.25:2]);
title('mode(infusion)/mode(baseline) IMM');

%% COMPARISON
figure; hold on
shadederrbar([1:3],ones(3,1)*nanmean(scop),ones(3,1)*SEM(scop,1),'k'); hold on
b = tmp_base; plot([1.1;2.1;3.1].*ones(3,4), b', '.b', 'MarkerSize', 20);
errorbar([1 2 3],nanmean(b,1),SEM(b,1),'b');
b = tmp; plot([1.1;2.1;3.1].*ones(3,4), b', '.r', 'MarkerSize', 20);
errorbar([1 2 3],nanmean(b,1),SEM(b,1),'r');

xlim([0.5 3.5]); xticks([1:3]); xticklabels({'aCSF','d1d2','glu'});
ylabel('voltage'); ylim([0 0.15]); yticks([0:0.05:0.15]);
title('mode(baseline) && mode(infusion) immobility only');

%%
figure; hold on
shadederrbar([1:3],ones(3,1)*nanmean(scop),ones(3,1)*SEM(scop,1),'k'); hold on
b = tmp; plot([1.1;2.1;3.1].*ones(3,4), b', '.r', 'MarkerSize', 20);
errorbar([1 2 3],nanmean(b,1),SEM(b,1),'r');

xlim([0.5 3.5]); xticks([1:3]); xticklabels({'aCSF','d1d2','glu'});
ylabel('voltage'); ylim([0 0.15]); yticks([0:0.05:0.15]);
title('mode(infusion) vs mode(mAChR) immobility only');