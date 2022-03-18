%% PLOT FFT w/o subtraction + DOWNSAMPLE
fig = figure; fig.Position(3) = 1375;
subplot(1,3,1); hold on
ds = 50;
plot([flog(f == 0.5) flog(f == 0.5)],[0 1],'--k'); 
plot([flog(f == 4) flog(f == 4)],[0 1],'--k'); 
shadederrbar(flog(1:ds:end), nanmean(norm_ach((1:ds:end),:),2), SEM(norm_ach((1:ds:end),:),2), [0.05 0.75 0.45]);
shadederrbar(flog(1:ds:end), nanmean(norm_da((1:ds:end),:),2), SEM(norm_da((1:ds:end),:),2), 'm');
shadederrbar(flog(1:ds:end), nanmean(norm_gfp((1:ds:end),:),2), SEM(norm_gfp((1:ds:end),:),2), 'k');
shadederrbar(flog(1:ds:end), nanmean(norm_tdt((1:ds:end),:),2), SEM(norm_tdt((1:ds:end),:),2), 'r');
legend({'0.5Hz','4Hz','ACh','','DA','','GFP','','tdt'})
xlabel('Frequency'); 
xlim([-1 flog(f == 50)]); xticks([-2:2]); xticklabels({'0.01','0.1','1','10','100'});
ylabel('Power (a.u.)');
axis('square');

% PLOT SUBTRACTION + SHADEDERRBAR
sub_ach = []; for x = 1:size(norm_ach,2); sub_ach(:,x) = norm_ach(:,x) - nanmean(norm_gfp,2); end % Subtract avg FFT for GFP
sub_da = []; for x = 1:size(norm_da,2); sub_da(:,x) = norm_da(:,x) - nanmean(norm_tdt,2); end % Subtract avg FFT for tdTomato

subplot(1,3,2); hold on
plot([-2 2],[0 0],'--k');
plot([flog(f == 0.5) flog(f == 0.5)],[-0.2 0.3],'--k'); 
plot([flog(f == 4) flog(f == 4)],[-0.2 0.3],'--k'); 
shadederrbar(flog(1:ds:end), nanmean(sub_ach((1:ds:end),:),2), SEM(sub_ach((1:ds:end),[1:4]),2), [0.05 0.75 0.45]); hold on
shadederrbar(flog(1:ds:end), nanmean(sub_da((1:ds:end),[5:8]),2), SEM(sub_da((1:ds:end),[5:8]),2), 'm'); hold on
legend({'line','0.5 Hz','4 Hz','ACh','DA'})
xlabel('Frequency'); 
xlim([-1 flog(f == 50)]); xticks([-2:2]); xticklabels({'0.01','0.1','1','10','100'});
title('FFT subtracting GFP'); ylabel('relative power (a.u.)');
axis('square');

% AUC post-subtraction
auc = [];
f_sub = 10.^(flog); % Regenerate frequency vector from log(freq)
r_14 = [find(f_sub == 0.5):find(f_sub == 4)]; % AUC from [1 4] Hz
for x = 1:size(sub_ach,2)
    auc(x,1) = trapz(sub_ach(r_14,x))/length(r_14);
    auc(x,2) = trapz(sub_da(r_14,x))/length(r_14);
end
%
subplot(1,3,3); hold on
errorbar(nanmean(auc)',SEM(auc,1)','.k','MarkerSize',20);
plot([1.15; 1.85].*ones(2,10),auc','.m','MarkerSize',20);
plot([1.15; 1.85].*ones(2,10),auc',':k','MarkerSize',20);
xlim([0.5 2.5]); xticks([1:2]); xticklabels({'ACh','DA'});
ylabel('relative power (a.u.)'); ylim([0 0.5]);
title(sprintf('ACh DA (n = %d mice)',size(auc,1)));
axis('square');