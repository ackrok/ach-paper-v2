%% adapted from getPhaseMat
%% Extract data
good_rew = [10:16,20,22:26,28:29,32:35,40:42,44:46]; % good rew
beh = modAChDA(good_rew);
beh = beh([5 7 13:19 23:25]);
[align_full, t_rew] = plot_fp2event(beh, [0 2], 0);

%% (1) Filter the signal into the frequency range of the oscillation
% Butterworth filter, set to filter a band of frequencies around the oscillation
passband = [0.1 5]; % in BZ case, theta was 8-10Hz in experiment so filter from 6-13Hz
Fs = 50;
filterFreq = passband/(Fs*0.5); 
[b, a] = butter(3, filterFreq, 'bandpass');

%%
fpMat = [];
phaseMat = [];
trialID = [];
trialRec = [];
win = 2; % seconds around event

for x = 1:length(beh); fpID = 2; 
%% (2) Convert the signal to phases
% Hilbert function in Matlab, then the angle function. 
% This will give you the signal in phases from -pi to pi.
signal = beh(x).FP{fpID} - nanmean(beh(x).FP{fpID}); % baseline adjust 
phases = angle(hilbert(filtfilt(b, a, signal))); % calculate phase of lfp

% make matrices for lfp and phases
rew = beh(x).reward./Fs;
rew(isnan(align_full{x,fpID}(1,:))) = []; % remove unrewarded trials

immOn = beh(x).onRest./Fs; immOff = beh(x).offRest./Fs; % immobility periods
immMid = ((immOff - immOn)/2) + immOn; % midpoint of immobility period

movOn = beh(x).on./Fs; % movement onset

events = [immMid(:);movOn(:);rew(:)];

sessionMat = nan(length(events), Fs*(2*win)+1);
phaseSessionMat = nan(length(events), Fs*(2*win)+1);

%% (3) Organize the signal into trials, centered on the reward time
% You should end up with a nTrial x nTimeBin matrix of phases after this step.
% As a side note, it's important to do steps 1-2 before step 3, to avoid edge
% effects from filtering smaller time segments of data.
trials = 1:length(events);
count = 0;
for iii = 1:length(events)
    count = count + 1;
    nTrial = trials(iii);
    index = round(Fs*events(iii));
    % add data from 3 second before and 3 seconds after
    sessionMat(nTrial,:) = signal((index - Fs*win):(index + Fs*win));
    phaseSessionMat(count,:) = phases((index - Fs*win):(index + Fs*win));
end

%% Store data
fpMat = [fpMat; sessionMat];
phaseMat = [phaseMat; phaseSessionMat];
trialID = [trialID; [1*ones(length(immMid),1); 2*ones(length(movOn),1); 3*ones(length(rew),1)]];
trialRec = [trialRec; x*ones(length(trials), 1)];

end

fprintf('phase reset analysis done! \n');
%% adapted from runPhaseResetTest
achda_phaseResetTest
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
