%% adj processFP 210319
%
% CHANGES:
%   - subtracting lowest 5% of dF/F values after processing
% ORDER
%   - after baseline/crunch/filter/downsample
%
% fpfinal = {};
for x = 1:length(raw)
    for y = 1:size(raw(x).fp,2)
        FPfinal = raw(x).fp(:,y);
        prcval = prctile(FPfinal, 1);
        FPfinal = FPfinal - prcval; % Subtract lower 5% of dF/F values
        % fpfinal{x,y} = FPfinal;
        raw(x).fp_prc1(:,y) = FPfinal;
    end
end

%% adj processFP 210318
%
% CHANGES: 
%   - added "crunch" adjustment to normalize pk2pk noise from 0 to 1
%       divide dF/F values by mean pk2pk noise across each sweep
%   - change fitType = 'interp', to interpolate points instead of fitting poly1
%
% ORDER:
%   (1) baselineFP: dF/F from RAW photometry signal (Voltage > dF/F)
%   (2) crunchFP to normalize noise levels to be between 0 and 1 dF/F
%   (3) filter <10Hz, downsample to 50Hz

params = struct;
params.lpCut = 20; params.filtOrder = 8;
params.rawFs = 2000; params.dsRate = 40; params.dsType = 2;
params.interpType = 'linear'; params.fitType = 'interp'; params.basePrc = 5; params.winSize = 10; params.winOv = 0;
h = waitbar(0,'processFP');
for x = 1:length(raw)
    for y = 1:size(raw(x).rawfp,2)
        rawFP = raw(x).rawfp(:,y);
        
%         bin = 0.01; % in seconds
%         tmp = [];
%         for aa = 1:[(length(rawFP)/rawFs)/bin]
%             bin_fp = rawFP([aa*(bin*rawFs)-(bin*rawFs-1)]:[aa*(bin*rawFs)]); % rawFP segment of length = bin
%             tmp(aa,1) = max(bin_fp); % maximum of bin segment
%             tmp(aa,2) = min(bin_fp); % minimun of bin segment
%         end
%         maxmin = mean(tmp(:,1) - tmp(:,2)); % pk2pk for each bin
%         % pk2pk_max = mean(tmp(:,1)); pk2pk_min = mean(tmp(:,2)); % mean pk-min, pk-max across all bins
%         raw(x).noiseV(y) = maxmin;
          % FPcrunch = (FPbase - pk2pk_min) ./ (maxmin/adjPk);  % min max normalization: noise is from 0 to 1
          
        [FPbase,baseline] = baselineFP(rawFP,params.interpType,params.fitType,params.basePrc,params.winSize,params.winOv,params.rawFs);
        % FPbase = (FPbase./100).*baseline; % alternatively, FPbase = raw_fp - baseline
        FPfilt = filterFP(FPbase,params.rawFs,params.lpCut,params.filtOrder,'lowpass');
        FPfinal = downsampleTLab(FPfilt,params.dsRate,params.dsType); Fs = params.rawFs/params.dsRate;
        raw(x).fp(:,y) = FPfinal;
        prcval = prctile(FPfinal, 1);
        FPfinal = FPfinal - prcval; % Subtract lower 5% of dF/F values
        raw(x).fp_prc1(:,y) = FPfinal;
        raw(x).params = params;
    end
    waitbar(x/length(raw),h);
end
close(h); fprintf('Done! \n');

%% adj processFP 210315
%
% CHANGES:
%   - adjusted parameters: basePrc = 5; winSize = 10; winOv = 0;
%
% ORDER:
%   (1) filterFP, (2) downsample, (3) baselineFP
%
% lpCut = 10; filtOrder = 8;
% rawFs = 2000; dsRate = 40; dsType = 2;
% interpType = 'linear'; fitType = 'interp'; basePrc = 5; winSize = 10; winOv = 0;
% 
% params.FP.lpCut = lpCut; params.FP.filtOrder = filtOrder;
% params.dsRate = dsRate; params.dsType = dsType; 
% params.FP.interpType = interpType; params.FP.fitType = fitType;
% params.FP.basePrc = basePrc; params.FP.winSize = winSize; params.FP.winOv = winOv;
% raw(1).params = params;
% 
% fpfinal = {};
% for x = 1:length(raw)
%     for y = 1:size(raw(x).rawfp,2)
%         rawFP = raw(x).rawfp(:,y);
%         nbFP = filterFP(rawFP,rawFs,lpCut,filtOrder,'lowpass');
%         nbFP = downsampleTLab(nbFP,dsRate,dsType); Fs = rawFs/dsRate;
%         [FP,baseline] = baselineFP(nbFP,interpType,fitType,basePrc,winSize,winOv,Fs);
%         %raw(x).fp(:,y) = FP;
%         fpfinal{x,y} = FP;
%     end
% end
% fprintf('Done!\n');

%%
figure;
for x = 1:length(raw)
    for y = 1:size(raw(x).fp,2)
        sp((3*(x-1))+y) = subplot(length(raw), 3, (3*(x-1))+y);
        plot(raw(x).fp(:,y)); 
        title(sprintf('%s-%d',raw(x).rec,y));
    end
end; linkaxes(sp,'y');

%%
Fs = 50;            % Sampling frequency                    
T = 1/Fs;             % Sampling period      
Y = fft(raw(x).fp(:,2));
L = length(Y);        % Length of signal
t = (0:L-1)*T;        % Time vector
f = (0:length(Y)-1)*Fs/length(Y);

plot(f,abs(Y))
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Magnitude')

%%
lpCut = 10; filtOrder = 8;
rawFs = 2000; dsRate = 40; dsType = 2;
interpType = 'linear'; fitType = 'line'; basePrc = 10; winSize = 15; winOv = 10;

basePrc = [1 5 10 50 75];
figure;
for aa = 1:length(basePrc)
x = 1; y = 3;
        rawFP = raw(x).rawfp(:,y);
        nbFP = filterFP(rawFP,rawFs,lpCut,filtOrder,'lowpass');
        %nbFP = filterFP(nbFP,rawFs,0.5,filtOrder,'highpass');
        nbFP = downsampleTLab(nbFP,dsRate,dsType); Fs = rawFs/dsRate;
        [FP,baseline] = baselineFP(nbFP,interpType,fitType,basePrc(aa),winSize,winOv,Fs);

        sp(aa) = subplot(length(basePrc),1,aa);
        plot(FP); title(sprintf('basePrc = %d | winSize = %d | winOv = %d',basePrc(aa),winSize,winOv)); grid on
end; linkaxes(sp,'y');

%% PLOT
figure
for x = 1:length(raw)
    for y = 1:size(raw(x).fp,2)
        sp((3*(x-1))+y) = subplot(length(raw),3,(3*(x-1))+y);
        plot(raw(x).fp(:,y)); ylim([0 1]);
    end
end
linkaxes(sp,'y');
