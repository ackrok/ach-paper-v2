%% new_ach new_da new_scop s raw

%% extract during IMMOBILITY + during INFUSION window
%glu all AK190_210909: new(x).onRest(z) < 6000000
%aCSFm d1d2 dhbe use 9000000: new(x).onRest(z) < 9000000
for x = 1:length(new)
    tmp = [];
    for z = 1:length(new(x).onRest) % Iterate over each immobility period
        % if new(x).onRest(z) < 6000000; continue; end % Extract only from t = +30min post infusion start time
        % if length(tmp) > 3000000; continue; end
        % tmp = [tmp; new(x).demod(new(x).onRest(z):new(x).offRest(z))]; % Extract and concatenate immobility periods
        tmp = [tmp; new(x).rawFP{2}(new(x).onRest(z):new(x).offRest(z))]; % Extract and concatenate immobility periods
%         for z = 1:length(new(x).on)
%         tmp = [tmp; new(x).rawFP{1}(new(x).on(z):new(x).off(z))]; % Extract and concatenate immobility periods
    end
    new(x).da_imm = tmp;
end

%% FFT
%new = new_ach([1:8 13 14 11 12]);
p1_mat = [];
for x = 1:length(new)
    vec = [new(x).da_imm]; 
    Fs = new(x).rawFs;
    
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
    p1_mat(:,x) = [movmean(P1,500)];
end
fprintf('FFT done! \n');

%% Normalize FFT
tmp = [];
r = [find(f == 0.1):find(f == 100)]; % Restrict to [0.01 100]
flog = log10(f(r));
f_sub = 10.^(flog); % Regenerate frequency vector from log(freq)
for x = 1:size(p1_mat,2)
    a = log10(p1_mat(r,x));
    vec_norm = (a - a(end))./(a(1) - a(end)); 
    % vec_norm = normalize(log10(p1_mat(r,x)),'range'); % Normalize range from [0.01 100], scaling so range covers [0 1]
    % vec_norm = normalize(log10(p1_mat(r,x)),'zscore'); % Normalize z-score
    % vec_norm = normalize(log10(p1_mat(r,x)),'range');
    % vec_norm = normalize(vec_norm,'zscore');
    tmp(:,x) = vec_norm;
end
norm_d1d2 = tmp;
fprintf('Normalization done! \n');

%% PLOT FFT w/o subtraction
figure; hold on
plot([flog(f == 0.5) flog(f == 0.5)],[0 1],'--k'); 
plot([flog(f == 4) flog(f == 4)],[0 1],'--k'); 
plot(flog, nanmean(norm(:,[1:4]),2), 'Color', [0.05 0.75 0.45]);
plot(flog, nanmean(norm(:,[5:8]),2), 'm');
plot(flog, nanmean(norm(:,[9:12]),2), 'Color', [0.85 0.35 0.1]);
plot(flog, nanmean(norm(:,[13:15]),2), 'b');
plot(flog, nanmean(norm_gfp,2), 'k');
legend({'0.5Hz','4Hz','aCSF','d1d2','glu','dhbe','GFP'})
xlabel('Frequency'); 
xlim([-2 2]); xticks([-2:2]); xticklabels({'0.01','0.1','1','10','100'});
% xlim([-1 1]); xticks([-1:1]); xticklabels({'0.1','1','10'});

%% PLOT FFT w/o subtraction + DOWNSAMPLE
figure; hold on
ds = 50;
plot([flog(f == 0.5) flog(f == 0.5)],[0 1],'--k'); 
plot([flog(f == 4) flog(f == 4)],[0 1],'--k'); 
shadederrbar(flog(1:ds:end), nanmean(norm((1:ds:end),[1:4]),2), SEM(norm((1:ds:end),[1:4]),2), [0.05 0.75 0.45]);
shadederrbar(flog(1:ds:end), nanmean(norm((1:ds:end),[5:8]),2), SEM(norm((1:ds:end),[5:8]),2), 'm');
shadederrbar(flog(1:ds:end), nanmean(norm((1:ds:end),[9:12]),2), SEM(norm((1:ds:end),[9:12]),2), [0.85 0.35 0.1]);
shadederrbar(flog(1:ds:end), nanmean(norm((1:ds:end),[13:15]),2), SEM(norm((1:ds:end),[13:15]),2), 'b');
shadederrbar(flog(1:ds:end), nanmean(norm_gfp((1:ds:end),:),2), SEM(norm_gfp((1:ds:end),:),2), 'k');
legend({'0.5Hz','4Hz','aCSF','','d1d2','','glu','','dhbe','','GFP',''})
xlabel('Frequency'); 
xlim([-1 flog(f == 50)]); xticks([-2:2]); xticklabels({'0.01','0.1','1','10','100'});
ylabel('Power (a.u.)');
axis('square');

%% PLOT SUBTRACTION
sub_gfp = []; for x = 1:size(norm,2); sub_gfp(:,x) = norm(:,x) - nanmean(norm_gfp,2); end % Subtract avg FFT for mAChR antagonist
figure; hold on
plot([-2 2],[0 0],'--k');
plot([0 0],[-0.2 0.3],'--k'); plot([0.6021 0.6021],[-0.2 0.3],'--k'); 
plot(flog, nanmean(sub_gfp(:,[1:4]),2), 'Color', [0.05 0.75 0.45]);
plot(flog, nanmean(sub_gfp(:,[5:8]),2), 'Color', 'm');
plot(flog, nanmean(sub_gfp(:,[9:12]),2), 'Color', [0.85 0.35 0.1]);
plot(flog, nanmean(sub_gfp(:,[13:15]),2), 'Color', 'b');
legend({'line','0.5 Hz','4 Hz','aCSF','d1d2','glu'})
xlabel('Frequency'); xlim([-2 2]); xticks([-2:2]); xticklabels({'0.01','0.1','1','10','100'});
title('FFT subtracting GFP');

%% PLOT SUBTRACTION + SHADEDERRBAR
sub_gfp = []; for x = 1:size(norm,2); sub_gfp(:,x) = norm(:,x) - nanmean(norm_gfp,2); end % Subtract avg FFT for mAChR antagonist
%%
ds = 50;
figure; hold on
plot([-2 2],[0 0],'--k');
plot([flog(f == 0.5) flog(f == 0.5)],[-0.2 0.3],'--k'); 
plot([flog(f == 4) flog(f == 4)],[-0.2 0.3],'--k'); 
shadederrbar(flog(1:ds:end), nanmean(sub_gfp((1:ds:end),[1:4]),2), SEM(sub_gfp((1:ds:end),[1:4]),2), [0.05 0.75 0.45]); hold on
shadederrbar(flog(1:ds:end), nanmean(sub_gfp((1:ds:end),[5:8]),2), SEM(sub_gfp((1:ds:end),[5:8]),2), 'm'); hold on
shadederrbar(flog(1:ds:end), nanmean(sub_gfp((1:ds:end),[9:12]),2), SEM(sub_gfp((1:ds:end),[9:12]),2), [0.85 0.35 0.1]); hold on
shadederrbar(flog(1:ds:end), nanmean(sub_gfp((1:ds:end),[13:15]),2), SEM(sub_gfp((1:ds:end),[13:15]),2), 'b'); hold on
legend({'line','0.5 Hz','4 Hz','aCSF','','d1d2','','glu','','nAChR'})
xlabel('Frequency'); 
xlim([-1 flog(f == 50)]); xticks([-2:2]); xticklabels({'0.01','0.1','1','10','100'});
title('FFT subtracting GFP'); ylabel('Power relative to GFP (a.u.)');

%% AUC post-subtraction
auc = [];
f_sub = 10.^(flog); % Regenerate frequency vector from log(freq)
r_14 = [find(f_sub == 0.5):find(f_sub == 4)]; % AUC from [1 4] Hz
for x = 1:size(sub_gfp,2)
    auc(x) = trapz(sub_gfp(r_14,x))/length(r_14);
end
auc = [auc, nan];
auc = reshape(auc, [4 4]);
nAn = 4;

% figure; violinplot(auc);  xticklabels({'aCSF','d1d2','glu'})
% hold on; plot([0.5 3.5],[0 0],'--k');
% title('AUC post-subtraction');
% AUC PLOT
a = auc; a = abs(a);
% a = auc./nanmean(auc(:,1)); % normalize to aCSF

figure; hold on
% plot([0.5 4.5],[1 1],'--k');
errorbar(nanmean(a)',SEM(a,1)','.k','MarkerSize',20);
plot([1.15; 1.85; 2.85; 3.85].*ones(4,4),a','.k','MarkerSize',20);
plot([1.15; 1.85].*ones(2,nAn),a(:,[1 2])',':m','MarkerSize',20);
plot([1.15; 2.85].*ones(2,nAn),a(:,[1 3])',':r','MarkerSize',20);
plot([1.15; 3.85].*ones(2,nAn),a(:,[1 4])',':b','MarkerSize',20);
xlim([0.5 4.5]); xticks([1:4]); xticklabels({'aCSF','D1R/D2R','NMDA/AMPA','DHbE'});

p = [];
[~,p(1)] = ttest(a(:,1));
[~,p(2)] = ttest(a(:,1),a(:,2));
[~,p(3)] = ttest(a(:,1),a(:,3));
[~,p(4)] = ttest(a(:,1),a(:,4));
title(sprintf('aCSF/d1d2: %1.3f | aCSF/glu: %1.3f | aCSF/DHbE: %1.3f',p(2),p(3),p(4)))

%% AUC PLOT normalize to aCSF
a = auc./auc(:,1);
