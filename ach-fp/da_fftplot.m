%fix norm_tdt
for x = 1:size(norm_gfp,2)
    a = norm_gfp(:,x);
    norm_gfp(:,x) = (a - a(125000))./(a(1) - a(125000));
end

% norm = [norm_da, norm_d2r(:,[3,5:8]), norm_lesion(:,[1 2 3 5 6])];
%% PLOT FFT w/o subtraction + DOWNSAMPLE
fig = figure; fig.Position(3) = 1375; 

subplot(1,3,1); hold on
ds = 50;
plot([flog(f == 0.5) flog(f == 0.5)],[0 1],'--k'); 
plot([flog(f == 4) flog(f == 4)],[0 1],'--k'); 
% shadederrbar(flog(1:ds:end), nanmean(norm_da_dual((1:ds:end),[1 3 5 9]),2), SEM(norm_da_dual((1:ds:end),[1 3 5 9]),2), 'm');
% shadederrbar(flog(1:ds:end), nanmean(norm_da_dual((1:ds:end),[2:2:12]),2), SEM(norm_da_dual((1:ds:end),[2:2:12]),2), 'g');
% shadederrbar(flog(1:ds:end), nanmean(norm_da_dms((1:ds:end),:),2), SEM(norm_da_dual((1:ds:end),:),2), 'b');
% shadederrbar(flog(1:ds:end), nanmean(norm_tdt((1:ds:end),:),2), SEM(norm_tdt((1:ds:end),:),2), 'k');
% legend({'0.5Hz','4Hz','DA DLS','','DA DMS','','DMS','','tdTomato'})

shadederrbar(flog(1:ds:end), nanmean(norm_da((1:ds:end),[1:4]),2), SEM(norm_da((1:ds:end),[1:4]),2), 'k');
shadederrbar(flog(1:ds:end), nanmean(norm_da((1:ds:end),[5:8]),2), SEM(norm_da((1:ds:end),[5:8]),2), 'm');
legend({'0.5Hz','4Hz','DA aCSF','','antag'})

xlabel('Frequency'); 
xlim([-1 flog(f == 50)]); xticks([-2:2]); xticklabels({'0.01','0.1','1','10','100'});
ylabel('Power (a.u.)'); ylim([0 1]); yticks([0:0.2:1]);
axis('square');

% PLOT SUBTRACTION + SHADEDERRBAR
% norm = [norm_da_dual, norm_da_dms];
% sub_gfp = []; for x = 1:size(norm,2); sub_gfp(:,x) = norm(:,x) - nanmean(norm_gfp(:,:),2); end % Subtract avg FFT for mAChR antagonist
sub_gfp = []; for x = 1:size(norm_da,2); sub_gfp(:,x) = norm_da(:,x) - nanmean(norm_da(:,[8]),2); end % Subtract avg FFT for mAChR antagonist
%
ds = 50;
subplot(1,3,2); hold on
plot([-2 2],[0 0],'--k');
plot([flog(f == 0.5) flog(f == 0.5)],[-0.2 0.3],'--k'); 
plot([flog(f == 4) flog(f == 4)],[-0.2 0.3],'--k'); 
% shadederrbar(flog(1:ds:end), nanmean(sub_gfp((1:ds:end),[1 3 5 9]),2), SEM(sub_gfp((1:ds:end),[1 3 5 9]),2), 'm'); 
% shadederrbar(flog(1:ds:end), nanmean(sub_gfp((1:ds:end),[2:2:12]),2), SEM(sub_gfp((1:ds:end),[2:2:12]),2), 'g'); 
% shadederrbar(flog(1:ds:end), nanmean(sub_gfp((1:ds:end),[13:19]),2), SEM(sub_gfp((1:ds:end),[13:19]),2), 'b'); 
% legend({'','0.5Hz','4Hz','DA DLS','','DA DMS','','DMS'})
shadederrbar(flog(1:ds:end), nanmean(sub_gfp((1:ds:end),[1:4]),2), SEM(sub_gfp((1:ds:end),[1:4]),2), 'k'); 
shadederrbar(flog(1:ds:end), nanmean(sub_gfp((1:ds:end),[5:8]),2), SEM(sub_gfp((1:ds:end),[5:8]),2), 'm'); 
legend({'','0.5Hz','4Hz','DA aCSF','','antag'})

xlabel('Frequency'); 
xlim([-1 flog(f == 50)]); xticks([-2:2]); xticklabels({'0.01','0.1','1','10','100'});
title('FFT subtracting tdTomato'); ylabel('Power relative to tdTomato (a.u.)');
axis('square');

% AUC post-subtraction
auc = [];
f_sub = 10.^(flog); % Regenerate frequency vector from log(freq)
r_14 = [find(f_sub == 0.5):find(f_sub == 4)]; % AUC from [1 4] Hz
for x = 1:size(sub_gfp,2)
    auc(x) = trapz(sub_gfp(r_14,x))/length(r_14);
end
auc(auc < 0) = 0;
% a = [auc(1:2:11)', auc(2:2:12)'];
a = [auc(1:4); auc(5:8)]';

subplot(1,3,3); hold on
%violinplot(auc_c);
plot(a', '--.k', 'MarkerSize', 20); xlim([0.5 2.5]); xticks([1 2]);
errorbar([0.75 2.22],nanmean(a), SEM(a,1), '.k', 'MarkerSize', 20);
% errorbar(cellfun(@nanmean, auc_c), cellfun(@nanstd, auc_c)./sqrt(cellfun(@length, auc_c)), '.k', 'MarkerSize', 20);
xticklabels({'DA control','d2r antag','lesion'})
hold on; plot([0.5 3.5],[0 0],'--k');
title('AUC post-subtraction'); axis('square');
ylim([-0.1 0.5]);