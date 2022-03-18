full = struct;
x = 1; full(x).FPnames = 'ACh';    full(x).beh = modAChDA; full(x).beh([17:25 27 47]) = [];
x = 2; full(x).FPnames = 'antag';  full(x).beh = raw_ant;
x = 3; full(x).FPnames = 'KO';     full(x).beh = raw_ko;   full(x).beh([1 2 4 6 11]) = [];
x = 4; full(x).FPnames = 'GFP';    full(x).beh = raw_gfp; 

%% Extract ACh dF/F values for IMMOBILITY
for s = 1:length(full)
    beh = full(s).beh;
    mat = struct;
    for x = 1:length(beh)
        mat(x).rec = beh(x).rec;
        mat(x).FPnames = beh(x).FPnames;
        mat(x).fp_imm = []; mat(x).fp_mov = []; % Initialize output vector

        if iscell(beh(x).onRest); y_end = length(beh(x).onRest); % Generalized for single-sweep and multi-sweep recordings
        else; y_end = 1; end
        for y = 1:y_end

            if iscell(beh(x).onRest)
                fp = beh(x).fp(:,y); % Extract photometry signal for this sweep
                idx_imm = extractEventST([1:length(fp)]',beh(x).onRest{y},beh(x).offRest{y},1); % Index of samples during immobility
                idx_mov = extractEventST([1:length(fp)]',beh(x).on{y},beh(x).off{y},1); % Index of samples during locomotion
            else
                if isfield(beh,'FP'); fp = beh(x).FP{1}; % Extract photometry signal
                elseif isfield(beh,'fp'); fp = beh(x).fp(:,y); end
                idx_imm = extractEventST([1:length(fp)]',beh(x).onRest,beh(x).offRest,1); % Index of samples during immobility
                idx_mov = extractEventST([1:length(fp)]',beh(x).on,beh(x).off,1); % Index of samples during locomotion
                if isfield(beh,'reward'); if ~isempty(beh(x).reward)
                    idx_rew = extractEventST([1:length(fp)]', floor(beh(x).reward), floor(beh(x).reward)+beh(x).Fs, 1);
                    idx_imm = idx_imm(~ismember(idx_imm,idx_rew)); % Remove indices of samples during reward window
                    idx_mov = idx_mov(~ismember(idx_mov,idx_rew)); % Remove indices of samples during reward window
                    end; end
            end
            % fp = fp - nanmean(fp);
            mat(x).fp_imm = [mat(x).fp_imm; fp(idx_imm)];
            mat(x).fp_mov = [mat(x).fp_mov; fp(idx_mov)];
        end
    end
    full(s).fp_imm = {mat.fp_imm}; full(s).fp_mov = {mat.fp_mov};
end
fprintf('Done \n');

% COMPUTE STATS
for s = 1:length(full)
    mu = []; med = []; sig = []; sig_mov = []; cv = [];
    for x = 1:length(full(s).fp_imm)
        mu(x) = nanmean(full(s).fp_imm{x});
        % med(x) = median(full(s).fp_imm{x});
        % sig(x) = SEM(mat(x).fp_imm,1);
        % sig(x) = nanstd(full(s).fp_imm{x});
        % sig_mov(x) = nanstd(full(s).fp_mov{x});
        cv(x,1) = nanstd(full(s).fp_imm{x})./nanmean(full(s).fp_imm{x});
        cv(x,2) = nanstd(full(s).fp_mov{x})./nanmean(full(s).fp_mov{x});
    end
    full(s).mu = mu; 
    % full(s).med = med; full(s).sig = sig; full(s).sig_mov = sig_mov;
    full(s).cv = cv;
end

% ADJUST by animal
cv_imm = cell(1, length(full)); cv_mov = cell(1,length(full));
for s = 1:length(full)
    tmp = {}; for x = 1:length(full(s).beh); tmp{x} = strtok(full(s).beh(x).rec,'-'); end
    uni = unique(tmp);
    % sig_an = []; sig_mov_an = [];
    cv_an = [];
    for x = 1:length(uni); idx = strcmp(tmp, uni{x});
        % sig_an(x) = nanmean(full(s).sig(idx));
        % sig_mov_an(x) = nanmean(full(s).sig_mov(idx));
        cv_an(x,:) = nanmean(full(s).cv(idx,:),1);
    end
    % full(s).sig_an = sig_an; full(s).sig_mov_an = sig_mov_an;
    full(s).cv_an = cv_an;
    cv_imm{s} = full(s).cv_an(:,1); cv_mov{s} = full(s).cv_an(:,2); 
end
fprintf('Done \n');

%%
% plot_me = cv_imm;
plot_me = {cv_imm{1}, cv_mov{1}};
mu = cellfun(@nanmean,plot_me);
std = cellfun(@nanstd,plot_me)./sqrt(cellfun(@length,plot_me));

figure; hold on
violinplot(plot_me);
errorbar(mu, std, 'k');
xlim([0.5 4.5]); xticks([1:4]); xticklabels({full.FPnames});
% p = []; for x = 2:4; [~,p(x)] = ttest2(plot_me{1},plot_me{x}); end
p = []; for x = 2:4; p(x) = ranksum(plot_me{1},plot_me{x}); end
title(sprintf('CV of IMM: ach/ant %1.3f || ach/ko %1.3f || ach/gfp %1.3f',p(2),p(3),p(4)));
% ylabel('Standard Deviation of ACh IMM Distribution'); ylim([0 4])
ylabel('CV (coefficient of variation)'); ylim([0 1]); yticks([0:0.25:1]);
axis('square');

%% TEST mu, std 
% mu = []; med = []; sig = [];
% for x = 1:length(full(s).fp_imm)
%     mu(x) = nanmean(full(s).fp_imm{x});
%     med(x) = median(full(s).fp_imm{x});
%     % sig(x) = SEM(mat(x).fp_imm,1);
%     sig(x) = nanstd(full(s).fp_imm{x});
% end
% figure; hold on
% plot(med, '.b', 'MarkerSize', 40); errorbar(mu, sig, '.r'); ylim([-5 5])
% 
%% PLOT distributions for each recording
% figure;
% for x = 1:length(mat)
%     sp(x) = subplot(7,7,x);
%     histogram(mat(x).fp_imm,'BinWidth',0.2,'Normalization','probability');
%     title(sprintf('%s',beh(x).rec));
% end
% linkaxes(sp,'x');

%% ADJUST by animal: IMMOBILITY vs LOCOMOTION
s = 1;
tmp = {}; for x = 1:length(full(s).beh); tmp{x} = strtok(full(s).beh(x).rec,'-'); end
uni = unique(tmp);
sig_an = []; sig_mov_an = [];
for x = 1:length(uni); idx = strcmp(tmp, uni{x});
    sig_an(x) = nanmean(full(s).sig(idx));
    sig_mov_an(x) = nanmean(full(s).sig_mov(idx));
end
sig_an_mat = [sig_an(:), sig_mov_an(:)]; 

a = sig_an_mat; nAn = size(a,1);
figure; hold on
errorbar(nanmean(a)',SEM(a,1)','.k','MarkerSize',20);
plot([1.25; 1.75].*ones(2,nAn),a','.:k','MarkerSize',20);
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'imm','mov'});
ylabel('standard deviation (s)'); ylim([0 7]); yticks([0:7]);
axis('square');
[~,p] = ttest(a(:,1),a(:,2));
title(sprintf('STD - I:%1.2f L:%1.2f \n p-value: I/L - %1.3f',nanmean(a(:,1)),nanmean(a(:,2)),p)); 