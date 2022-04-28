%% adapted from immPhaseResetTest
%% (4) Run the statistical test of choice on the columns of the matrix from step 3
% This will tell you if the phase distribution is uniform or not.
% Circular Statistics toolbox from Philipp Berens -- circ_rtest. 
% This is the Raleigh test, which is the strongest test for unimodal deviations from uniformity. 
% The documentation does a really good job of explaining circular statistics and how to use each function.

%% (5) Show the distributions of phase in time
% For this, make a 2d histogram of the data in phase and time (using histcounts2 in matlab), 
% and look for the distribution of phases in time at the reward time. 
% There should be a visible peak in phases if there is a phase reset. 
% Also, if the distribution happens to be bimodal, not unimodal, then it 
% is necessary to run a different statistical test, such as the Omnibus test. 
% This is also in the circular statistics package

%% ALL STATISTICS
switch fpID
    case 1; lbl = {'ACh immobility','ACh onset','ACh reward'};
    case 2; lbl = {'DA immobility','DA onset','DA reward'};
end
id = {[find(trialID == 1)], [find(trialID == 2)], [find(trialID == 3)]};
colors = [1 0 0; 0 1 0; 0 0 1]; % b, r, g
N = 3;

circStd = cell(N,1); % circular standard deviation
circMean = cell(N,3); % mean direction, with upper and lower 95% CI
skew = cell(N,1); % skewness
kurt = cell(N,1); % excess kurtosis
p_ral = cell(N,1); % Rayleigh Test of Uniformity
p_omn = cell(N,1); % Omnibus Test of Uniformity
rVecLen = cell(N,1); % resultant vector mean length for the phases
for x = 1:N
    circStd{x} = circ_std(phaseMat(id{x},:)); % circular standard deviation for circular data 
    [circMean{x,1}, circMean{x,2}, circMean{x,3}] = circ_mean(phaseMat(id{x},:)); % mean direction for circular data
    skew{x} = circ_skewness(phaseMat(id{x},:)); % skewness
    [~, kurt{x}] = circ_kurtosis(phaseMat(id{x}, :)); % excess Kurtosis
    p_ral{x} = nan(1, size(phaseMat, 2));
    p_omn{x} = nan(1, size(phaseMat, 2));
    rVecLen{x} = nan(size(phaseMat, 2), 1); 
    for ii = 1:size(phaseMat,2)
        p_ral{x}(ii) = circ_rtest(phaseMat(id{x},ii)); % Rayleigh Test of Uniformity
        p_omn{x}(ii) = circ_otest(phaseMat(id{x},ii)); % Omnibus Test of Uniformity
        rVecLen{x}(ii) = circ_r(phaseMat(id{x}, ii)); % calculate resultant vector mean for the phases
    end
end

% %% While Rayleigh test is  p < 0.05, the mean phase and 95% confidence interval
dphi = 0.1;
phaseBins = linspace(-pi, pi, 64);
nPhiBins = length(phaseBins)-1;

dt = 0.015;
timeBins = -win:dt:win;

sd = 0.7;
alpha = 0.05; % threshold: while Rayleigh test is p < 0.05, find mean phase and 95% confidence intervals
time = linspace(-win, win, size(fpMat,2));
% note: win is set in achda_phaseReset1

% 95% CI for mean phase direction in circMean{:,[2 3]}
% significance of mean:
p_mean = nan(1, length(circMean{3,1}));
for ii = 1:length(circMean{3,1})
    p_mean(ii) = circ_mtest(phaseMat(id{3}, ii), circMean{3,1}(ii));
end

% %% Distribution of phase in time
% For this, make a 2d histogram of the data in phase and time (using histcounts2 in matlab), 
% and look for the distribution of phases in time at the reward time. 
% There should be a visible peak in phases if there is a phase reset. 
p_phase = cell(N,1);
for x = 1:N
    phase_match = phaseMat(id{x},:);
    p_phase{x} = histcounts2(reshape(phase_match', 1, []), ...
        repmat(time, 1, size(phase_match, 1)), phaseBins, timeBins, ...
        'Normalization', 'probability');
end

clc
fprintf('Done phase reset statistics! \n');

%% PLOT MEAN + STDEV
% figure;
% subplot(2,1,1); hold on
% for x = 1:N; plot(circMean{x,1}, 'color', colors(x,:), 'lineWidth', 2); end
% xlabel('Time (s)'); xticks((1:Fs/2:size(fpMat,2))); xticklabels(-win:0.5:win);
% ylabel('Mean'); ylim([-pi pi])
% title('phase');
% 
% subplot(2,1,2); hold on
% for x = 1:N; plot(circStd{x}, 'color', colors(x,:), 'lineWidth', 2); end
% xlabel('Time (s)'); xticks((1:Fs/2:size(fpMat, 2))); xticklabels(-win:0.5:win);
% ylabel('StDev');
% legend(lbl,'Location','southeast')

%% PLOT MEAN PHASE for each BEH
% figure; clearvars sp_mean
% for x = 1:N
%     sp_mean(x) = subplot(3,1,x); hold on
%     plot(circMean{x,1}, 'color', colors(x,:), 'lineWidth', 2); % mean direction
%     plot(circMean{x,2}, 'k', 'lineWidth', 0.5); % upper 95% confidence interval
%     plot(circMean{x,3}, 'k', 'lineWidth', 0.5); % lower 95% confidence interval
%     % xlabel('Time (s)'); 
%     xticks((1:Fs/2:size(fpMat,2))); xticklabels(-win:0.5:win); xlim([1 size(fpMat,2)]);
%     ylabel('Mean');
%     title(lbl{x});
% end
% linkaxes(sp_mean,'x'); linkaxes(sp_mean,'y');

%% Skewness, Kurtosis
% % Note that the skewness and kurtosis are not very useful, and are heavily
% % influenced by sample size. 
% figure;
% subplot(2,1,1); hold on
% for x = 1:N; plot(skew{x}, 'color', colors(x,:), 'lineWidth', 2); end
% xlabel('Time (s)'); 
% xticks((1:Fs/2:size(fpMat,2))); xticklabels(-win:0.5:win); xlim([1 size(fpMat,2)]);
% ylabel('Skewness'); title('Skewness');
% ylim([-1 1])
% 
% subplot(2,1,2); hold on
% for x = 1:N; plot(kurt{x}, 'color', colors(x,:), 'lineWidth', 2); end
% xlabel('Time (s)'); 
% xticks((1:Fs/2:size(fpMat,2))); xticklabels(-win:0.5:win); xlim([1 size(fpMat,2)]);
% ylabel('Excess Kurtosis'); title('Excess Kurtosis');
% legend(lbl)
% % ylim([-1 1])

%% Use the Rayleigh test (for unimodal) and omnibus test (unimodal, biomodal, and multimodal) 
% % to test for a uniform distribution
% figure;
% sp_test = [];
% sp_test(1) = subplot(2,1,1); hold on
% for x = 1:N; plot(p_ral{x}, 'color', colors(x,:), 'lineWidth', 2); end
% plot(1:length(p_ral{3}), 0.05*ones(1, length(p_ral{3})), 'color', [0.5 0.5 0.5])
% xlabel('Time (s)'); 
% xticks((1:Fs/2:size(fpMat,2))); xticklabels(-win:0.5:win); xlim([1 size(fpMat,2)]);
% ylabel('p'); ylim([0 1])
% legend(lbl)
% title('Rayleigh Test of Uniformity')
% 
% sp_test(2) = subplot(2,1,2); hold on
% for x = 1:N; plot(p_omn{x}, 'color', colors(x,:), 'lineWidth', 2); end
% plot(1:length(p_omn{3}), 0.05*ones(1, length(p_omn{3})), 'color', [0.5 0.5 0.5])
% xlabel('Time (s)'); 
% xticks((1:Fs/2:size(fpMat,2))); xticklabels(-win:0.5:win); xlim([1 size(fpMat,2)]);
% ylabel('p'); ylim([0 1])
% title('Omnibus Test of Uniformity')
% linkaxes(sp_test,'x');

%% DISTRIBUTION of phase in time
% fig = figure; % fig.Position = [20 20 1000 750];
% sp_dist = [];
% for x = 1:N
%     sp_dist(x) = subplot(3, 1, x); hold on
%     imagesc(imgaussfilt(p_phase{x}, sd));
%     ttVals = linspace(1, 166, length(circMean{x,1}));
%     rtest_indices = p_ral{x} < alpha & ~isnan(circMean{x,2});
%     
%     scatter(ttVals(rtest_indices), (circMean{x,1}(rtest_indices)+ pi)*9.87, 5, 'r', 'filled') % mean direction
%     scatter(ttVals(rtest_indices), (circMean{x,2}(rtest_indices) + pi)*9.87, 3, 'k', 'filled') % upper 95% confidence interval
%     scatter(ttVals(rtest_indices), (circMean{x,3}(rtest_indices) + pi)*9.87, 3, 'k', 'filled') % lower 95% confidence interval
%     
%     ylabel('phase'); yticks([1 floor(nPhiBins/2) nPhiBins]); yticklabels({'-\pi', '0', '\pi'})
%     xticks(linspace(1, length(timeBins)-1, 5)); xticklabels(-win:win)
%     xlim([1 length(timeBins)-1]);
%     colorbar; % clim([5e-5 12e-5])
%     title(lbl{x});
% end
% linkaxes(sp_dist,'x');

%% FULL FIGURE
fig = figure; fig.Position = [20 20 1300 750];
sp = [];
for x = 1:N
    sp(x,1) = subplot(5, 3, x);
    imagesc(imgaussfilt(p_phase{x}, sd))
    xlabel('Time'); 
    xticks(linspace(1, length(timeBins)-1, 5)); xticklabels(-win:win)
    xlim([1 size(p_phase{x}, 2)])
    ylabel('phase'); 
    yticks([1 floor(nPhiBins/2) nPhiBins]); yticklabels({'-\pi', '0', '\pi'})
    axis xy
    colorbar; % clim([0 1.5e-4]) % clim([9e-5 10.5e-5])
    title(lbl{x})
    
    sp(x,2) = subplot(5, 3, x+3);
    plot(circMean{x,1}, 'color', colors(x,:), 'lineWidth', 2)
    xlabel('Time (s)'); 
    xticks((1:Fs:size(fpMat,2))); xticklabels(-win:win)
    xlim([1 length(circMean{x,1})])
    ylabel('mean'); ylim([-4 4]);
    
    sp(x,3) = subplot(5, 3, x+6);
    plot(circStd{x}, 'color', colors(x,:), 'lineWidth', 2)
    xlabel('Time (s)'); 
    xticks((1:Fs:size(fpMat,2))); xticklabels(-win:win)
    xlim([1 length(circStd{x})])
    ylabel('\sigma'); ylim([0.5 1.5]);

    sp(x,4) = subplot(5, 3, x+9);
    plot(p_ral{x}, 'color', colors(x,:), 'lineWidth', 2)
    xlabel('Time (s)'); 
    xticks((1:Fs:size(fpMat,2))); xticklabels(-win:win)
    xlim([1 length(p_ral{x})])
    ylabel('p(uniform)'); ylim([0 1]);
    
    sp(x,5) = subplot(5, 3, x+12);
    plot(rVecLen{x}, 'color', colors(x,:), 'lineWidth', 2)
    xlabel('Time (s)'); 
    xticks((1:Fs:size(fpMat,2))); xticklabels(-win:win)
    xlim([1 length(rVecLen{x})])
    ylabel('Resultant Vector Length'); ylim([0 1]);
end
linkaxes(sp(:,1),'x'); 
linkaxes(sp(:,2),'y'); linkaxes(sp(:,3),'y'); linkaxes(sp(:,4),'y'); linkaxes(sp(:,5),'y');

