%% DA FFT
fig = figure; fig.Position(3) = 1375;
ds = 50;
subplot(1,3,1); hold on
shadederrbar(flog(1:ds:end), nanmean(norm((1:ds:end),[1:4]),2), SEM(norm((1:ds:end),[1:4]),2), 'm'); hold on
shadederrbar(flog(1:ds:end), nanmean(ant(1:ds:end,:),2), SEM(ant(1:ds:end,:),2), 'k');
xlabel('Frequency'); xlim([-2 2]); xticks([-2:2]); xticklabels({'0.01','0.1','1','10','100'});
title('FFT - DOPAMINE'); ylabel('Power (a.u.)'); axis('square');

% FFT subtraction
% sub = []; for x = 1:4
%     sub(:,x) = norm((1:ds:end),x) - nanmean(ant(1:ds:end,:),2); end
subplot(1,3,2); hold on
plot([-2 2],[0 0],'--k');
plot([log10(0.5) log10(0.5)],[-0.2 0.3],'--k'); 
plot([log10(4) log10(4)],[-0.2 0.3],'--k'); 
shadederrbar(flog(1:ds:end), nanmean(sub,2), SEM(sub,2), 'm');
legend({'line','0.5 Hz','4 Hz','aCSF-d1d2'})
xlabel('Frequency'); ylabel('Power relative to antag (a.u.)');
xlim([-2 2]); xticks([-2:2]); xticklabels({'0.01','0.1','1','10','100'});
title('FFT DA subtracting DA-antag'); axis('square');

% AUC post-subtraction
auc = [];
f_sub = 10.^(flog); % Regenerate frequency vector from log(freq)
r_14 = [find(f_sub == 0.5):find(f_sub == 4)]; % AUC from [1 4] Hz
for x = 1:size(tmp,2)
    auc(x) = trapz(tmp(r_14,x))/length(r_14);
end
auc = [auc, nan];
auc = reshape(auc, [4 3]);
a = auc(:,[1 3]);


subplot(1,3,3); hold on
plot([0.5 2.5],[0 0],'--k');
errorbar(nanmean(a)',SEM(a,1)','.k','MarkerSize',20);
plot([1.15; 1.85].*ones(2,4),a','--m','MarkerSize',20);
xlim([0.5 2.5]); xticks([1:2]); xticklabels({'aCSF','DHbE'});
ylim([-0.1 0.3]); yticks([-0.1:0.1:0.3]); ylabel('Power relative to antag (a.u.)');
[~,p] = ttest(a(:,1),a(:,2));
title(sprintf('acsf vs dhbe p = %1.3f',p)); axis('square');

%%
figure; hold on
plot([0.5 3.5],[0 0],'--k');
errorbar(nanmean(a)',SEM(a,1)','.k','MarkerSize',20);
plot([1.15; 1.85].*ones(2,4),a','.--m','MarkerSize',20);
errorbar(3, nanmean(auc), SEM(auc,2)','.r','MarkerSize',20);
plot([2.85].*ones(1,6),auc','.r','MarkerSize',20);
xlim([0.5 3.5]); xticks([1:3]); xticklabels({'aCSF','DHbE','b2flox'});
ylim([-0.1 0.3]); yticks([-0.1:0.1:0.3]); ylabel('Power relative to antag (a.u.)');
[~,p] = ttest(a(:,1),a(:,2)); [~,p2] = ttest2(a(:,1),auc);
title(sprintf('acsf vs dhbe p = %1.3f || vs beta2 p = %1.3f',p,p2)); axis('square');