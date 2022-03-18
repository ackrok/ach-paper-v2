
%% extract during IMMOBILITY + during INFUSION window
for x = 1:length(new)
    tmp = [];
    for z = 1:length(new(x).onRest)
        tmp = [tmp; new(x).rawFP{1}(new(x).onRest(z):new(x).offRest(z))]; % Extract and concatenate immobility periods
    end
    new(x).imm = tmp;
end
%%
% raw_all = struct;
for x = 1:length(raw)
    n = 1+length(raw_all);
    raw_all(n).rec = raw(x).rec;
    raw_all(n).FPnames = raw(x).FPnames{1};
    rawfp_imm = []; rawfp_mov = []; 
    for y = 1:size(raw(x).rawfp,2)
        for z = 1:length(raw(x).onRest{y})
            rawfp_imm = [rawfp_imm; raw(x).rawfp(raw(x).onRest{y}(z)*(2000/50):raw(x).offRest{y}(z)*(2000/50),y)]; 
        end
        for z = 1:length(raw(x).on{y})
            rawfp_mov = [rawfp_mov; raw(x).rawfp(raw(x).on{y}(z)*(2000/50):raw(x).off{y}(z)*(2000/50),y)]; 
        end
    end
    raw_all(n).rawfp_imm = rawfp_imm;
    raw_all(n).rawfp_mov = rawfp_mov;
end
%%
p1_mat = [];
for x = 1:length(raw)
    vec = [raw(x).rawfp_imm]; 
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
    p1_mat(:,x) = [movmean(P1,500)];
end
fprintf('FFT done! \n');

%% Normalize FFT
tmp = [];
r = [find(f == 0.1):find(f == 100)]; % Restrict to [0.01 100]
flog = log10(f(r));
f_sub = 10.^(flog); % Regenerate frequency vector from log(freq)
for x = 1:size(p1_mat,2)
    vec_norm = normalize(log10(p1_mat(r,x)),'range'); % Normalize range from [0.01 100], scaling so range covers [0 1]
    tmp(:,x) = vec_norm;
end
norm = tmp;

%% adjust by ANIMAL
norm = [norm_ach, norm_ko, norm_ant, norm_gfp];
tmp = {}; for x = 1:length(raw_imm); tmp{x} = strtok(raw_imm(x).rec,'-'); end
uni = unique(tmp);
norm_an = []; for x = 1:length(uni); idx = find(strcmp(tmp,uni{x}));
    norm_an(:,x) = nanmean(norm(:,idx),2); end

norm_ach_an = [norm_an(:,[5:9]),norm_ach(:,[23:32])];
norm_ko_an = norm_an(:,[10:14]);
norm_gfp_an = norm_an(:,[15:17]);

%% Plot NORMALIZED + SHADEDERRBAR
figure; hold on
ds = 50;
plot([flog(f == 0.5) flog(f == 0.5)],[0 1],'--k'); 
plot([flog(f == 4) flog(f == 4)],[0 1],'--k'); 
% shadederrbar(flog(1:ds:end), nanmean(norm_ach((1:ds:end),:),2), SEM(norm_ach((1:ds:end),:),2), [0.05 0.75 0.45]);
shadederrbar(flog(1:ds:end), nanmean(norm_ach_an((1:ds:end),:),2), SEM(norm_ach_an((1:ds:end),:),2), [0.05 0.75 0.45]);
shadederrbar(flog(1:ds:end), nanmean(norm_ko_an((1:ds:end),:),2), SEM(norm_ko_an((1:ds:end),:),2), 'r');
shadederrbar(flog(1:ds:end), nanmean(norm_ant((1:ds:end),:),2), SEM(norm_ant((1:ds:end),:),2), 'b');
shadederrbar(flog(1:ds:end), nanmean(norm_gfp_an((1:ds:end),:),2), SEM(norm_gfp_an((1:ds:end),:),2), 'k');
%legend({'0.5Hz','4Hz','ACh','','KO','','ANT','','GFP',''})
xlabel('Frequency'); 
xlim([-1 log10(50)]); xticks([-2:2]); xticklabels({'0.01','0.1','1','10','100'});
ylabel('Power (a.u.)');
title('ACh n = 15 || KO n = 5 || ANT n = 5 || GFP n = 3'); axis('square');

%% PLOT SUBTRACTION + SHADEDERRBAR
norm = [norm_ach_an, norm_ko_an, norm_ant];
sub_gfp = []; for x = 1:size(norm,2); sub_gfp(:,x) = norm(:,x) - nanmean(norm_gfp_an,2); end % Subtract avg FFT 

figure; hold on
ds = 50;
plot([-2 2],[0 0],'--k');
plot([flog(f == 0.5) flog(f == 0.5)],[-0.2 0.3],'--k'); 
plot([flog(f == 4) flog(f == 4)],[-0.2 0.3],'--k'); 
% shadederrbar(flog(1:ds:end), nanmean(sub_gfp((1:ds:end),[1:22]),2), SEM(sub_gfp((1:ds:end),[1:22]),2), [0.05 0.75 0.45]); hold on
% shadederrbar(flog(1:ds:end), nanmean(sub_gfp((1:ds:end),[23:40]),2), SEM(sub_gfp((1:ds:end),[23:40]),2), 'r'); hold on
% shadederrbar(flog(1:ds:end), nanmean(sub_gfp((1:ds:end),[41:45]),2), SEM(sub_gfp((1:ds:end),[41:45]),2), 'b'); hold on
shadederrbar(flog(1:ds:end), nanmean(sub_gfp((1:ds:end),[1:15]),2), SEM(sub_gfp((1:ds:end),[1:15]),2), [0.05 0.75 0.45]); hold on
shadederrbar(flog(1:ds:end), nanmean(sub_gfp((1:ds:end),[16:20]),2), SEM(sub_gfp((1:ds:end),[16:20]),2), 'r'); hold on
shadederrbar(flog(1:ds:end), nanmean(sub_gfp((1:ds:end),[21:25]),2), SEM(sub_gfp((1:ds:end),[21:25]),2), 'b'); hold on
legend({'line','0.5 Hz','4 Hz','ACh','','KO','','antag',''})
xlabel('Frequency'); xlim([-1 1]); xticks([-2:2]); xticklabels({'0.01','0.1','1','10','100'});
title('FFT subtracting GFP'); axis('square');
ylabel('Power relative to GFP (a.u.)'); ylim([-0.2 0.5]); yticks([-0.2:0.1:0.5]);

%% AUC post-subtraction, N = #mice
auc = [];
f_sub = 10.^(flog); % Regenerate frequency vector from log(freq)
r_14 = [find(f_sub == 0.5):find(f_sub == 4)]; % AUC from [1 4] Hz
for x = 1:size(sub_gfp,2)
    auc(x) = trapz(sub_gfp(r_14,x))/length(r_14);
end
% auc_cell = {auc(1:15),auc(16:20),auc(21:25)};
auc_cell = {auc(1:15)./nanmean(auc(1:15)),auc(16:20)./nanmean(auc(1:15)),auc(21:25)./nanmean(auc(1:15))};

figure; violinplot(auc_cell_norm);  
xticklabels({'ACh','KO','antag'})
hold on; plot([0.5 3.5],[0 0],'--k');
errorbar([1:3],cellfun(@nanmean, auc_cell),cellfun(@nanstd, auc_cell)./sqrt(cellfun(@length, auc_cell)),'k');
p = []; for x = 2:3; [~,p(x)] = ttest2(auc_cell{1},auc_cell{x}); end
title(sprintf('AUC post-subtract GFP \n RS p-value: wt/KO %1.4f | wt/ANT %1.4f',ranksum(auc_cell{1},auc_cell{2}),ranksum(auc_cell{1},auc_cell{3})));

%% AUC post-subtraction N = #mice
% tmp = {}; for x = 1:45; tmp{x} = strtok(raw_imm(x).rec,'-'); end
% uni = unique(tmp);
% auc_an = [];
% for x = 1:length(uni); idx = find(strcmp(tmp,uni{x}));
%     auc_an(x) = nanmean(auc(idx)); end
% auc_cell_an = {auc_an(1:5),auc_an(6:8),auc_an(9:13)};
% 
% figure; violinplot(auc_cell_an);  xticklabels({'ACh','KO','antag'})
% hold on; plot([0.5 3.5],[0 0],'--k');
% errorbar([1:3],cellfun(@nanmean, auc_cell_an),cellfun(@nanstd, auc_cell_an)./sqrt(cellfun(@length, auc_cell_an)),'k');
% p = []; for x = 1:3; [~,p(x)] = ttest(auc_cell_an{x}); end
% title(sprintf('AUC [0.5 4Hz] post-subtract GFP \n p-value: ACh %1.3f | KO %1.3f | ANT %1.3f',p(1),p(2),p(3)));
% ylim([-0.2 0.5]); yticks([-0.2:0.1:0.5]);

% a = auc;
% figure; hold on
% plot([0.5 3.5],[0 0],'--k');
% errorbar(nanmean(a)',SEM(a,1)','.k','MarkerSize',20);
% plot([1.15; 1.85; 2.85].*ones(3,4),a','.m','MarkerSize',20);
% xlim([0.5 3.5]); xticks([1:3]); xticklabels({'aCSF','D1R/D2R','NMDA/AMPA'});