%%
tmp = [];
h = waitbar(0, 'FFT + normalize');
for x = 1:length(raw)
    vec = [raw(x).imm]; 
    Fs = raw(x).rawFs;
    
    needL = 2500*Fs;
    vec = repmat(vec,[ceil(needL/length(vec)) 1]);
    vec = vec(1:needL);
    T = 1/Fs;               % Sampling period
    L = length(vec);        % Length of signal
    vec(isnan(vec)) = [];
    fftACh = fft(vec);      % Discrete Fourier Transform of photometry signal
    P2 = abs(fftACh/L);     % Two-sided spectrum P2
    P1 = P2(1:L/2+1);       % Single-sided spectrum P1 based on P2 and even-valued signal length L
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;     % Frequency domain vector
    P1 = medfilt1(P1);      % Median filter initial FFT
    P1 = movmean(P1,500);
    
    r = [find(f == 0.1):find(f == 100)]; % Restrict to [0.01 100]
    flog = log10(f(r));
    f_sub = 10.^(flog); % Regenerate frequency vector from log(freq)
    
    vec_norm = normalize(log10(P1(r)),'range'); % Normalize range from [0.01 100], scaling so range covers [0 1]
    tmp(:,x) = vec_norm;
    waitbar(x/length(raw),h);
end
fprintf('FFT done! \n'); close(h);

norm_dual = tmp;

%% adjust by ANIMAL
norm_ach_imm_an = []; 
norm_ach_mov_an = [];
tmp = {}; for x = 1:32; tmp{x} = strtok(raw_imm(x).rec,'-'); end
uni = unique(tmp);
for x = 1:length(uni); idx = find(strcmp(tmp,uni{x}));
    norm_ach_imm_an(:,x) = nanmean(norm_imm(:,idx),2); 
    norm_ach_mov_an(:,x) = nanmean(norm_mov(:,idx),2); 
end

norm_gfp_imm_an = []; 
norm_gfp_mov_an = [];
tmp = {}; for x = 56:67; tmp{x-55} = strtok(raw_imm(x).rec,'-'); end
uni = unique(tmp);
for x = 1:length(uni); idx = 55 + find(strcmp(tmp,uni{x}));
    norm_gfp_imm_an(:,x) = nanmean(norm_imm(:,idx),2); 
    norm_gfp_mov_an(:,x) = nanmean(norm_mov(:,idx),2); 
end

%% PLOT NORMALIZED
% fig = figure; fig.Position(3) = 1375;

figure; hold on
% subplot(1,3,1); hold on
ds = 50;
% plot([flog(f == 0.5) flog(f == 0.5)],[0 1],'--k'); 
% plot([flog(f == 4) flog(f == 4)],[0 1],'--k');
% tmp = [norm_ach_imm_an((1:ds:end),[10:12])];% , norm_dual((1:ds:end),[7 11])];
% shadederrbar(flog(1:ds:end), nanmean(tmp,2), SEM(tmp,2), [0.05 0.75 0.45]);
% shadederrbar(flog(1:ds:end), nanmean(norm_ach_imm_an((1:ds:end),:),2), SEM(norm_ach_imm_an((1:ds:end),:),2), [0.05 0.75 0.45]);
% shadederrbar(flog(1:ds:end), nanmean(norm_dms_an((1:ds:end),:),2), SEM(norm_dms_an((1:ds:end),:),2), 'c');
% shadederrbar(flog(1:ds:end), nanmean(norm_ach_mov_an((1:ds:end),:),2), SEM(norm_ach_mov_an((1:ds:end),:),2), 'b');
shadederrbar(flog(1:ds:end), nanmean(norm_gfp_imm_an((1:ds:end),:),2), SEM(norm_gfp_imm_an((1:ds:end),:),2), 'k');
% shadederrbar(flog(1:ds:end), nanmean(norm_gfp_mov_an((1:ds:end),:),2), SEM(norm_gfp_mov_an((1:ds:end),:),2), 'c');
shadederrbar(flog(1:ds:end), nanmean(norm_dual((1:ds:end),[1 3 5 9 11]),2), SEM(norm_ach_imm_an((1:ds:end),[1 3 5 9 11]),2), [0.05 0.75 0.45]);
shadederrbar(flog(1:ds:end), nanmean(norm_dual((1:ds:end),[2 4 6 8 12]),2), SEM(norm_ach_imm_an((1:ds:end),[2 4 6 8 12]),2), 'b');
% legend({'0.5Hz','4Hz','GFP','DLS','DMS'});
xlabel('Frequency'); 
xlim([-1 log10(50)]); xticks([-2:2]); xticklabels({'0.01','0.1','1','10','100'});
ylabel('Power (a.u.)');
title('ACh DLS n = 5 || DMS n = 5 (matched)'); axis('square');

%% PLOT SUBTRACTION
% norm = [norm_ach_imm_an, norm_dms_an];
norm = norm_dual([1 3 5 9 11, 2 4 6 8 12]);
sub_gfp = []; 
for x = 1:size(norm,2); sub_gfp(:,x) = norm(:,x) - nanmean(norm_gfp_imm_an,2); end % Subtract avg FFT 

figure; hold on
% subplot(1,3,2); hold on
ds = 50;
plot([-2 2],[0 0],'--k');
plot([flog(f == 0.5) flog(f == 0.5)],[-0.2 0.3],'--k'); 
plot([flog(f == 4) flog(f == 4)],[-0.2 0.3],'--k'); 
shadederrbar(flog(1:ds:end), nanmean(sub_gfp((1:ds:end),[1:5]),2), SEM(sub_gfp((1:ds:end),[1:5]),2), [0.05 0.75 0.45]);
shadederrbar(flog(1:ds:end), nanmean(sub_gfp((1:ds:end),[6:10]),2), SEM(sub_gfp((1:ds:end),[6:10]),2), 'c');
legend({'line','0.5 Hz','4 Hz','DLS','','DMS'})
xlabel('Frequency'); xlim([-1 log10(50)]); xticks([-2:2]); xticklabels({'0.01','0.1','1','10','100'});
title('FFT subtracting GFP'); axis('square');
ylabel('Power relative to GFP (a.u.)'); ylim([-0.2 0.5]); yticks([-0.2:0.1:0.5]);

%% AREA UNDER CURVE
auc = [];
f_sub = 10.^(flog); % Regenerate frequency vector from log(freq)
r_14 = [find(f_sub == 0.5):find(f_sub == 4)]; % AUC from [1 4] Hz
for x = 1:size(sub_gfp,2)
    auc(x,1) = trapz(sub_gfp(r_14,x))/length(r_14);
end

a = {auc(1:15),auc(16:20)};
figure; hold on
% subplot(1,3,3); hold on
violinplot(a);
errorbar(cellfun(@nanmean,a)',(cellfun(@nanstd,a)./sqrt(cellfun(@length,a)))','.k','MarkerSize',20);
ylabel('relative power (a.u.)'); ylim([0 0.5]); yticks([0:0.1:0.5]);
[~,p] = ttest2(a{1},a{2});
title(sprintf('ranksum %1.4f, ttest2 %1.4f',ranksum(a{1},a{2}),p));

% a = auc; nAn = size(sub_gfp,2);
% subplot(1,3,3); hold on
% errorbar(nanmean(a)',SEM(a,1)','.k','MarkerSize',20);
% plot([1.15; 1.85].*ones(2,nAn),a','.-g','MarkerSize',20);
% xlim([0.5 2.5]); xticks([1 2]); xticklabels({'ACh IMM','ACh MOV'});
% ylabel('relative power (a.u.)'); ylim([0 0.5]); yticks([0:0.1:0.5]);
% [~,p] = ttest(a(:,1),a(:,2));
% title(sprintf('signrank %1.4f, ttest %1.4f',signrank(a(:,1),a(:,2)),p));
% axis('square');
% movegui(fig,'center');