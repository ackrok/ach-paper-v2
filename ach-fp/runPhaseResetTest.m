% look at phase variance

% load data
region = 'ca1';
figDir = fullfile('D:', 'BuzsakiRinzel', 'Figures', 'LFP');
thisDir = pwd;

colors = cbrewer('qual', 'Set2', 8, 'cubic');

filename = [basename '.phaseMat.mat'];

load(filename)

%%  if want to look at all sessions
jump = condID < 7;
run = condID > 6;

%% Visualize LFP

figure
hold on
for ii = find(jump)'
    plot(LFPMat(ii,:) + ii*6000, 'k')
end
box off
plot(1875*ones(1,2), [0 ii*6000 + 2000], 'r')
xticks((1:625:size(LFPMat, 2))); xticklabels(-1.5:0.5:1)
title('Jump')

figure
hold on
for ii = find(run)'
    plot(LFPMat(ii,:) + ii*6000, 'k')
end
plot(1875*ones(1,2), [9.5e5 ii*6000 + 2000], 'r')
box off
xticks((1:625:size(LFPMat, 2))); xticklabels(-1.5:0.5:1)
title('Run')


%% Run general statistics

figure
subplot(2,1,1)
plot(circ_mean(phaseMat(jump,:)),'color', colors(1,:), 'lineWidth', 2)
hold on
plot(circ_mean(phaseMat(run,:)), 'color', colors(2,:),'lineWidth', 2)
box off
xticks((1:625:size(LFPMat, 2))); xticklabels(-1.5:0.5:1)
xlabel('Time (s)')
ylabel('Mean')
legend('Jump', 'Run')
ylim([-pi pi])

subplot(2,1,2)
[std_jump] = circ_std(phaseMat(jump, :));
[std_run] = circ_std(phaseMat(run, :));
plot(std_jump,'color', colors(1,:), 'lineWidth', 2)
hold on
plot(std_run, 'color', colors(2,:),'lineWidth', 2)
box off
xticks((1:625:size(LFPMat, 2))); xticklabels(-1.5:0.5:1)
xlabel('Time (s)')
ylabel('St. Dev')
% legend('Jump', 'Run')
% ylim([1.35 1.45])

%%
figure
subplot(2,1,1)
[mj, ulj, llj] = circ_mean(phaseMat(jump, :));
[mr, ulr, llr] = circ_mean(phaseMat(run, :));
plot(mj, 'color', colors(1,:), 'lineWidth', 2)
hold on
plot(ulj, 'k', 'lineWidth', 0.5)
plot(llj, 'k', 'lineWidth', 0.5)
xticks((1:625:size(LFPMat, 2))); xticklabels(-1.5:0.5:1)
xlabel('Time (s)')
ylabel('Mean')
subplot(2,1,2)
plot(mr, 'color', colors(2,:), 'lineWidth', 2)
hold on
plot(ulr, 'k', 'lineWidth', 0.5)
plot(llr, 'k', 'lineWidth', 0.5)
xticks((1:625:size(LFPMat, 2))); xticklabels(-1.5:0.5:1)
xlabel('Time (s)')
ylabel('Mean')

% Note that the skewness and kurtosis are not very useful, and are heavily
% influenced by sample size. 
%%

figure
subplot(2,1,1)
plot(circ_skewness(phaseMat(jump,:)),'color', colors(1,:), 'lineWidth', 2)
hold on
plot(circ_skewness(phaseMat(run,:)), 'color', colors(2,:),'lineWidth', 2)
box off
xticks((1:625:size(LFPMat, 2))); xticklabels(-1.5:0.5:1)
xlabel('Time (s)')
ylabel('Skewness')
legend('Jump', 'Run')
ylim([-1 1])

subplot(2,1,2)
[~, k_jump] = circ_kurtosis(phaseMat(jump, :));
[~, k_run] = circ_kurtosis(phaseMat(run, :));
plot(k_jump,'color', colors(1,:), 'lineWidth', 2)
hold on
plot(k_run, 'color', colors(2,:),'lineWidth', 2)
box off
xticks((1:625:size(LFPMat, 2))); xticklabels(-1.5:0.5:1)
xlabel('Time (s)')
ylabel('Excess Kurtosis')
legend('Jump', 'Run')
ylim([-1 1])


%% Use the Rayleigh test (for unimodal) and omnibus test (unimodal, biomodal, and multimodal) to test for a uniform distribution
p_jr = nan(1, size(phaseMat, 2));
p_rr = nan(1, size(phaseMat, 2));
for ii = 1:size(phaseMat, 2)
p_jr(ii) = circ_rtest(phaseMat(jump,ii));
p_rr(ii) = circ_rtest(phaseMat(run,ii));
end

p_jo = nan(1, size(phaseMat, 2));
p_ro = nan(1, size(phaseMat, 2));
for ii = 1:size(phaseMat, 2)
p_jo(ii) = circ_otest(phaseMat(jump,ii));
p_ro(ii) = circ_otest(phaseMat(run,ii));
end

figure
subplot(2,1,1)
plot(p_jr,'color', colors(1,:), 'lineWidth', 2)
hold on
plot(p_rr, 'color', colors(2,:),'lineWidth', 2)
plot(1:length(p_jr), 0.05*ones(1, length(p_jr)), 'color', [0.5 0.5 0.5])
box off
xticks((1:625:size(LFPMat, 2))); xticklabels(-1.5:0.5:1)
xlabel('Time (s)')
ylabel('p')
% legend('Jump', 'Run')
ylim([0 1])
title('Rayleigh Test of Uniformity')

subplot(2,1,2)
plot(p_jo,'color', colors(1,:), 'lineWidth', 2)
hold on
plot(p_ro, 'color', colors(2,:),'lineWidth', 2)
plot(1:length(p_jr), 0.05*ones(1, length(p_jr)), 'color', [0.5 0.5 0.5])
box off
xticks((1:625:size(LFPMat, 2))); xticklabels(-1.5:0.5:1)
xlabel('Time (s)')
ylabel('p')
% legend('Jump', 'Run')
ylim([0 1])
title('Omnibus Test of Uniformity')

%% Plot, while Rayleigh test is  p < 0.05, the mean phase and 95% confidence interval
dphi = 0.1;
% phaseBins = -pi:dphi:pi+ 0.1; % 
phaseBins = linspace(-pi, pi, 64);
nPhiBins = length(phaseBins)-1;

dt = 0.015;
timeBins = -1.5:dt:1;

sd = 0.7;
alpha = 0.05;
time = linspace(-1.5, 1, 3126);

jt = find(condID < 7);

phaseJump = phaseMat(jt, :);

p_phaseJ = histcounts2(reshape(phaseJump', 1, []), repmat(time, 1, size(phaseJump, 1)), phaseBins, timeBins, 'Normalization', 'probability');

% 95% confidence intervals
[mj, ulj, llj] = circ_mean(phaseMat(jump, :));
[mr, ulr, llr] = circ_mean(phaseMat(run, :));
% significance of mean
p_mean = nan(1, length(mj));
for ii = 1:length(mj)
    p_mean(ii) = circ_mtest(phaseMat(jump, ii), mj(ii));
end

fig = figure;
subplot(2, 1, 1)
imagesc(imgaussfilt(p_phaseJ, sd))
ttVals = linspace(1, 166, length(mj));
rtest_indices = p_jr < alpha & ~isnan(ulj);
hold on
scatter(ttVals(rtest_indices), (mj(rtest_indices)+ pi)*9.87, 5, 'r', 'filled')
hold on
scatter(ttVals(rtest_indices), (ulj(rtest_indices) + pi)*9.87, 3, 'k', 'filled')
scatter(ttVals(rtest_indices), (llj(rtest_indices) + pi)*9.87, 3, 'k', 'filled')
ylabel('phase')
xlabel('time')
xticks(linspace(1, length(timeBins)-1, 6)); xticklabels(-1.5:0.5:1.1)
yticks([1 floor(nPhiBins/2) nPhiBins]); yticklabels({'-\pi', '0', '\pi'})
colorbar
% clim([5e-5 12e-5])
box off
fig.Position = [20 20 1000 750];
title('Jump')

jt = find(condID > 6);
phaseJump = phaseMat(jt, :);

p_phaseR = histcounts2(reshape(phaseJump', 1, []), repmat(time, 1, size(phaseJump, 1)), phaseBins, timeBins, 'Normalization', 'probability');


subplot(2, 1, 2)
imagesc(imgaussfilt(p_phaseR, sd))
ttVals = linspace(1, 166, length(mr));
rtest_indices = p_rr < alpha & ~isnan(ulr);
hold on
scatter(ttVals(rtest_indices), (mr(rtest_indices)+ pi)*9.87, 5, 'r', 'filled')
hold on
scatter(ttVals(rtest_indices), (ulr(rtest_indices) + pi)*9.87, 3, 'k', 'filled')
scatter(ttVals(rtest_indices), (llr(rtest_indices) + pi)*9.87, 3, 'k', 'filled')
ylabel('phase')
xlabel('time')
xticks(linspace(1, length(timeBins)-1, 6)); xticklabels(-1.5:0.5:1.1)
yticks([1 floor(nPhiBins/2) nPhiBins]); yticklabels({'-\pi', '0', '\pi'})
colorbar
clim([5e-5 12e-5])
box off
fig.Position = [20 20 1000 750];
title('Run')

fig = figure;
subplot(5, 2, 1)
imagesc(imgaussfilt(p_phaseJ, sd))
ttVals = linspace(1, 166, length(mj));
rtest_indices = p_jr < alpha & ~isnan(ulj);
hold on
scatter(ttVals(rtest_indices), (mj(rtest_indices)+ pi)*9.87, 5, 'r', 'filled')
hold on
scatter(ttVals(rtest_indices), (ulj(rtest_indices) + pi)*9.87, 3, 'k', 'filled')
scatter(ttVals(rtest_indices), (llj(rtest_indices) + pi)*9.87, 3, 'k', 'filled')
ylabel('phase')
xlabel('time')
xticks(linspace(1, length(timeBins)-1, 6)); xticklabels(-1.5:0.5:1.1)
yticks([1 floor(nPhiBins/2) nPhiBins]); yticklabels({'-\pi', '0', '\pi'})
colorbar
clim([5e-5 12e-5])
box off
% fig.Position = [20 20 1000 750];
title('Jump')

jt = find(condID > 6);
phaseJump = phaseMat(jt, :);

p_phaseR = histcounts2(reshape(phaseJump', 1, []), repmat(time, 1, size(phaseJump, 1)), phaseBins, timeBins, 'Normalization', 'probability');


subplot(5, 2, 3)
imagesc(imgaussfilt(p_phaseR, sd))
ttVals = linspace(1, 166, length(mr));
rtest_indices = p_rr < alpha & ~isnan(ulr);
hold on
scatter(ttVals(rtest_indices), (mr(rtest_indices)+ pi)*9.87, 5, 'r', 'filled')
hold on
scatter(ttVals(rtest_indices), (ulr(rtest_indices) + pi)*9.87, 3, 'k', 'filled')
scatter(ttVals(rtest_indices), (llr(rtest_indices) + pi)*9.87, 3, 'k', 'filled')
ylabel('phase')
xlabel('time')
xticks(linspace(1, length(timeBins)-1, 6)); xticklabels(-1.5:0.5:1.1)
yticks([1 floor(nPhiBins/2) nPhiBins]); yticklabels({'-\pi', '0', '\pi'})
colorbar
clim([5e-5 12e-5])
box off
fig.Position = [20 20 1000 750];
title('Run')


%% calculate resultant vector mean for the phases
Rlen = nan(size(phaseMat, 2), 1);

for ii = 1:size(phaseMat, 2)
    Rlen(ii) = circ_r(phaseMat(jump, ii));
end


%% Okay, now make presentation figures

fig = figure;
subplot(5, 2, 1)
imagesc(imgaussfilt(p_phaseJ, sd))
ylabel('phase')
xlabel('Time')
axis xy
xticks(linspace(1, length(timeBins)-1, 6)); xticklabels(-1.5:0.5:1.1)
yticks([1 floor(nPhiBins/2) nPhiBins]); yticklabels({'-\pi', '0', '\pi'})
% colorbar
% clim([0 1.5e-4])
% clim([9e-5 10.5e-5])
xlim([1 size(p_phaseJ, 2)])
box off
% fig.Position = [20 20 1000 750];
title('Jump')

subplot(5, 2, 3)
plot(stdJumpVar, 'k', 'lineWidth', 2)
% plot(std_jump, 'k', 'lineWidth', 2)
box off
xticks((1:625:size(LFPMat, 2))); xticklabels(-1.5:0.5:1)
xlim([1 length(stdJumpVar)])
% ylim([1.37 1.43])
xlabel('Time')
ylabel('\sigma')

subplot(5, 2, 5)
plot(p_jr, 'k', 'lineWidth', 2)
box off
xticks((1:625:size(LFPMat, 2))); xticklabels(-1.5:0.5:1)
xlim([1 length(stdJumpVar)])
% ylim([1.25 1.41])
xlabel('Time')
ylabel('p(uniform)')

subplot(5, 2, 7)
imagesc(imgaussfilt(p_phaseJ, sd))
ylabel('phase')
xlabel('Time')
axis xy
xticks(linspace(1, length(timeBins)-1, 6)); xticklabels(-1.5:0.5:1.1)
yticks([1 floor(nPhiBins/2) nPhiBins]); yticklabels({'-\pi', '0', '\pi'})
colorbar
% clim([8e-5 11e-5])
% clim([0 1.5e-4])
% clim([7e-5 11e-5])
xlim([1 size(p_phaseJ, 2)])
box off
% fig.Position = [20 20 1000 750];
title('Jump')

fig.Position = [20 20 500 800];
% fig.PaperPosition = [.25 .25 10 8];

subplot(5, 2, 9)
plot(circ_std(phaseMat(jump, :)), 'k', 'lineWidth', 2)
box off
xticks((1:625:size(LFPMat, 2))); xticklabels(-1.5:0.5:1)
% xlim([1 length(stdJumpVar)])
% ylim([1.37 1.43])
xlabel('Time')
ylabel('\sigma')

subplot(5, 2, 2)
plot(Rlen, 'k', 'lineWidth', 2)
box off
xticks((1:625:size(LFPMat, 2))); xticklabels(-1.5:0.5:1)
% xlim([1 length(stdJumpVar)])
% ylim([1.37 1.43])
xlabel('Time')
ylabel('Resultant Vector Length')

