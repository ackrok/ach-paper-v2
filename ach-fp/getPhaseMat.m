% get phases of LFP only

clear all, close all, clc

% folder to save to:
saveDir = fullfile('D:', 'BuzsakiRinzel', 'Data');

%% Get session folder names to loop through
thisDir = pwd;

d = dir;

d = d(~ismember({d.name}, {'.', '..'}));

index = zeros(1, length(d));
for ii = 1:length(d)
    index(ii) = isfolder(d(ii).name);
end
d = d(logical(index));

%% set up matrices % parameters
LFPMat = [];
phaseMat = [];
condID = [];
sessionID = [];

region = 'hpc';

Fs = 50;
filterFreq = [0.1 5]/(Fs*0.5); % lfp.samplingRate = 1250
[b, a] = butter(3, filterFreq, 'bandpass');

%% Loop through sessions and save data
for ii = 1:length(d)
    
    try % not every session was processed or has an .lfp file, so code won't run. Instead, skip these folders
        
        cd(d(ii).name)
        fprintf(['Session ' d(ii).name '\n'])
        
        dirName = pwd;
        basename = bz_BasenameFromBasepath(dirName);
        
        
        % load behavior, spikes, and lfp files
        load([basename '.behavior.mat']);
        
        spikes = bz_GetSpikes('basepath', dirName);
        sessionInfo = bz_getSessionInfo;
        if any(strcmp(spikes.region, region))
            lfp = bz_GetLFP(sessionInfo.ca1);
        end
        
        % calculate phase of lfp
        phases = angle(hilbert(filtfilt(b, a, double(lfp.data))));
        
        % make matrices for lfp and phases
        sessionMat = nan(length(behavior.events.trials), lfp.samplingRate*6+1);
        phaseSessionMat = nan(length(behavior.events.trials), lfp.samplingRate*6+1);
        
        rate = lfp.samplingRate;
        
        %% get good trials
        
        trials = 1:length(behavior.events.trials);
        
        count = 0;
        for iii = 1:length(trials)
            count = count + 1;
            ntrial = trials(iii);
            if behavior.events.trialConditions(ntrial) < 7 % check if trial is a jump trial
                index = round(rate*behavior.events.jumpTime(ntrial, 1));
            else % if not a jump trial, center trial on halfway point
                index = round(rate*mean(behavior.events.trialIntervals(ntrial, :)));
            end
            if ~isnan(index) % add data from 3 second before jump to 3 seconds after.
                sessionMat(ntrial, :) = lfp.data((index-rate*3):(index+rate*3));
                phaseSessionMat(count, :) = phases((index-rate*3):(index+rate*3));
            end
        end
        
        % store data
        LFPMat = [LFPMat; sessionMat];
        phaseMat = [phaseMat; phaseSessionMat];
        condID = [condID; behavior.events.trialConditions(trials)'];
        sessionID = [sessionID; ii*ones(length(trials), 1)];
        
        fprintf([basename ' added \n'])
        
        
    catch e %e is an MException struct
        
        fprintf(1,'The identifier was:\n%s \n',e.identifier );
        fprintf(1,'There was an error! The message was:\n%s \n',e.message);
        
        fprintf([basename ' did not work \n'])
    end
    
    cd(thisDir) % return to main directory
end

%% Save data

filename = [region 'phaseMat.mat'];

save(filename, 'LFPMat', 'phaseMat', 'condID', 'sessionID')
