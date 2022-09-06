behLoaded = menu('beh loaded into workspace already?','yes','no');

switch behLoaded
    case 2
        %% Select .mat files you want to add to summary data structu
        fPath = 'R:\tritsn01labspace\'; 
        [fName,fPath] = uigetfile([fPath,'*.mat'],'MultiSelect','On');
        if ~iscell(fName); fName = {fName}; end

        %% Extract data
        beh = struct;
        for z = 1:length(fName) 
            load(fullfile(fPath,fName{z})); 
            [an,j2] = strtok(fName{z},'_'); day = strtok(j2,'_');
            x = 1+length(beh);
            beh(x).rec = [an,'-',day]; 
            beh(x).site = 'DLS';
            beh(x).task = 'wheel';

            %% Photometry
            beh(x).Fs = data.gen.Fs; % Sampling frequency, in Hz
            beh(x).time = data.final.time; % Time vector
            beh(x).FP = data.final.FP;
            beh(x).nbFP = data.final.nbFP; % Photometry signal(s)
            beh(x).FPnames = data.final.FPnames; % Names of photometry signal(s)

            %% Movement
            if isfield(data.final,'vel') % If movement data exists
                beh(x).vel = data.final.vel; % Velocity signal
                beh(x).on = data.final.mov.onsets; beh(x).off = data.final.mov.offsets;                 % Movement onset/offset times in sampling freq (data.gen.Fs), NOT in seconds
                beh(x).onRest = data.final.mov.onsetsRest; beh(x).offRest = data.final.mov.offsetsRest; % Rest onset/offset times in sampling freq (data.gen.Fs), NOT in seconds
            else % Open field
                try 
                    beh(x).cam = data.final.cam.on; 
                catch
                    signal = data.acq.opto.trace;
                    sigEdge = data.gen.params.FP.sigEdge; 
                    rawFs = data.gen.params.acqFs; dsRate = data.gen.params.dsRate;
                    if sigEdge ~= 0
                        signal = signal((sigEdge*rawFs)+1:end-(sigEdge*rawFs));
                    end
                    camOn = getPulseOnsetOffset (signal, 0.5);
                    camOn_Fs = round(camOn/dsRate);
                    beh(x).cam = camOn_Fs;
                end
            end

            %% Lick/Reward
            if isfield(data.acq,'rew')
                % data = processReward(data, data.gen.params);
                beh(x).task = 'reward';
                beh(x).reward = data.final.rew.onset;    % Reward delivery time in sampling freq (data.gen.Fs), NOT in seconds
                beh(x).lick = data.final.lick.onset;        % Lick times in sampling freq (data.gen.Fs), NOT in seconds
                beh(x).lickVec = data.final.lick.trace;        % Lick trace
            end

            %%
            fprintf('Extracted from %s\n',fName{z});
        end
        if isempty(beh(1).Fs); beh(1) = []; end
end
%%
if ~exist(beh)
    error('No variable called beh exists');
end

%% Variables
[plotOpt,ind] = listdlg('PromptString',{'Select All Plots to Generate',...
    'For Multiple: Hold Ctrl and Select'},'ListString',...
    {'DA/ACh cross-correlation',...
    'DA/ACh coherence/phase'});


%%
if ind == 0
    msgbox('Plotting Aborted');
else
    for plotOptNum = 1:length(plotOpt)
        choice = plotOpt(plotOptNum);
        switch choice
            case 1 % Cross-correlation
                %%
                [corr_achda, lags, shuff_achda, min_val, min_lag] = AK_corrFP(beh);
                nAn = size(corr_achda{1},2);
                %%
                fig = figure; fig.Position(3) = 1375;
                clr = {'r','g','b'};
                Fs = beh(1).Fs;
                subplot(1,3,1); hold on; 
                    plot([0 0],[-1 0.5],'--k');
                    shadederrbar(lags/Fs, nanmean(shuff_achda{1,2},2), nanmean(shuff_achda{1,1},2), 'k'); hold on
                    for z = 1:3
                        shadederrbar(lags/Fs, nanmean(corr_achda{z,1},2), SEM(corr_achda{z,1},2), clr{z});
                    end
                    xlabel('lag from DA (s)'); xlim([-1 1]);
                    ylabel('coefficient'); ylim([-1 0.5]); yticks([-1:0.25:1]);
                    title(sprintf('DA/ACh corr (n = %d mice)',nAn));
                    axis square
                %
                jit = [];
                for z = 1:3
                    j1 = z-0.25; j2 = z+0.25; jit(:,z) = j1 + (j2-j1).*rand(nAn,1); % jitter x-values for plotting
                end
                %
                subplot(1,3,2); hold on;
                    a = min_val;
                    plot(jit,a,'.k','MarkerSize',20); % plot raw values per animal
                    errorbar(nanmean(a,1), SEM(a,1),'.r', 'MarkerSize', 20);
                    xlim([0.5 3.5]); xticks([1:3]); xticklabels({'imm','mov','rew'}); 
                    ylabel('coefficient'); ylim([-1 0]); yticks([-1:0.2:1]);
                    p = anova1(a,[],'off');
                    title(sprintf('coeff: anova p = %1.2f',p));
                    axis square
                %    
                subplot(1,3,3); hold on;
                    min_lag(abs(min_lag) == 5) = nan; % if lag is +/-5s, then adjust to be nan
                    a = min_lag.*1000; % adjust lag at minimum to be in milliseconds, from seconds
                    plot(jit,a,'.k','MarkerSize',20); % plot raw values per animal
                    errorbar(nanmean(a,1), SEM(a,1),'.r', 'MarkerSize', 20);
                    xlim([0.5 3.5]); xticks([1:3]); xticklabels({'imm','mov','rew'}); 
                    ylabel('lag (ms)'); ylim([-250 0]); yticks([-250:50:0]);
                    p = anova1(a,[],'off');
                    title(sprintf('lag: anova p = %1.2f',p));
                    axis square
                %
                movegui(gcf,'center');
                
            case 2 % Coherence and phase
                %%
                [coher_achda, phase_achda, t, f, coher_shuff, phase_shuff] = AK_coherFP(beh);
                nAn = size(coher_achda{1},2);
                %%
                fig = figure; fig.Position([3 4]) = [1000 860];
                clr = {'r','g','b'};
                band = [0.5 4]; % frequency band: 0.5-4 Hz
                subplot(2,2,1); hold on
                    for z = 1:3 % iterate over behavioral states
                        % plot coherence magnitude, averaged across all
                        % animals for each behavioral state
                        shadederrbar(f, nanmean(coher_achda{z},2), SEM(coher_achda{z},2), clr{z});
                    end
                    plot([band(1) band(1)],[0 1],'--k'); plot([band(2) band(2)],[0 1],'--k'); % plot lines at 0.5, 4 Hz
                    shadederrbar(f, nanmean(coher_shuff{2},2), nanmean(coher_shuff{3},2) - nanmean(coher_shuff{2},2), 'k');
                    legend({'imm','','loc','','rew','','0.5Hz','4Hz'});
                    xlabel('frequency');
                    ylabel('coherence'); ylim([0 1]); yticks([0:0.2:1]);
                    title(sprintf('DA/ACh coherence (n = %d mice)',nAn));
                    axis square
                %
                subplot(2,2,2); hold on
                    for z = 1:3 % iterate over behavioral states
                        % plot phase (in degrees), averaged across all
                        % animals for each behavioral state
                        shadederrbar(f, nanmean(rad2deg(phase_achda{z}),2), SEM(rad2deg(phase_achda{z}),2), clr{z}); 
                    end
                    plot([band(1) band(1)],[-180 180],'--k'); plot([band(2) band(2)],[-180 180],'--k'); % plot lines at 0.5, 4 Hz
                    shadederrbar(f, nanmean(rad2deg(phase_shuff{2}),2), nanmean(rad2deg(phase_shuff{3}),2) - nanmean(rad2deg(phase_shuff{2}),2), 'k');
                    % legend({'imm','','loc','','rew','','0.5Hz','4Hz'});
                    xlabel('frequency');
                    ylabel('degrees'); ylim([-180 180]); yticks([-180:90:180]);
                    title('coherence phase'); 
                    axis square
                %
                r = [6:42]; % [~,r2] = min(abs(f - 2)); % range for 0.5-4Hz
                coher_avg = []; phase_avg = [];
                for z = 1:3
                    % coher_avg(:,z) = coher_achda{z}(r2,:); phase_avg(:,z) = phase_achda{z}(r2,:); % value at 2Hz
                    % coher_avg(:,z) = nanmean(coher_achda{z}(r,:)); phase_avg(:,z) = nanmean(phase_achda{z}(r,:)); % average within frequency band
                    coher_avg(:,z) = median(coher_achda{z}(r,:)); phase_avg(:,z) = median(phase_achda{z}(r,:)); % median within frequency band
                end
                jit = [];
                for z = 1:3
                    j1 = z-0.25; j2 = z+0.25; jit(:,z) = j1 + (j2-j1).*rand(nAn,1); % jitter x-values for plotting
                end
                %
                subplot(2,2,3); hold on
                    a = coher_avg;
                    plot(jit,a,'.k','MarkerSize',20); % plot raw values per animal
                    errorbar(nanmean(a,1), SEM(a,1),'.r', 'MarkerSize', 20);
                    xlim([0.5 3.5]); xticks([1:3]); xticklabels({'imm','mov','rew'}); 
                    ylabel('coherence'); ylim([0 1]); yticks([0:0.2:1]);
                    p = anova1(a,[],'off');
                    title(sprintf('coherence: anova p = %1.2f',p));
                    axis square
                %
                subplot(2,2,4); hold on
                    a = rad2deg(phase_avg);
                    plot(jit,a,'.k','MarkerSize',20); % plot raw values per animal
                    errorbar(nanmean(a,1), SEM(a,1),'.r', 'MarkerSize', 20);
                    xlim([0.5 3.5]); xticks([1:3]); xticklabels({'imm','mov','rew'}); 
                    ylabel('phase (deg)'); ylim([0 180]); yticks([0:45:180]);
                    p = anova1(a,[],'off');
                    title(sprintf('phase: anova p = %1.2f',p));
                    axis square
                %
                movegui(gcf,'center');
                % 
        end
    end
end