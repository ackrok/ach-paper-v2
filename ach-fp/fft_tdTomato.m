%%
fPath = 'C:\Users\Anya\Desktop\FP_LOCAL\'; 
[fName,fPath] = uigetfile([fPath,'*.mat'],'MultiSelect','On');
if ~iscell(fName); fName = {fName}; end
%%
new = struct;
h = waitbar(0, 'new processing');
for y = 1:length(fName)
    load(fullfile(fPath,fName{y})); % Load data file
    x = 1 + length(new);
    new(x).rec = dataRaw.mouserec; fprintf('%s ... ',new(x).rec);
    new(x).site = 'DLS';
    new(x).rx = 'tdTomato';
    %% Pull photometry data
    new(x).rawFP = dataRaw.acq.FP{2};
    new(x).rawFs = dataRaw.gen.acqFs;
    waitbar(y/length(fName),h);
    fprintf('DONE.\n');
end
close(h);

%%
fPath = 'C:\Users\Anya\Desktop\FP_LOCAL\'; 
[fName,fPath] = uigetfile([fPath,'*.mat'],'MultiSelect','On');
if ~iscell(fName); fName = {fName}; end
%%
h = waitbar(0, 'new processing');
for y = 1:length(fName)
    load(fullfile(fPath,fName{y})); % Load data file
    x = find(strcmp({new.rec},data.mouserec));
    dsRate = data.gen.acqFs/data.gen.Fs;
    new(x).onRest = data.final.mov.onsetsRest*dsRate;
    new(x).offRest = data.final.mov.offsetsRest*dsRate;
    waitbar(y/length(fName),h);
    fprintf('DONE.\n');
end
close(h);
%%
for x = 6:9
    tmp = [];
    h = waitbar(0, 'new processing');
    for z = 1:length(new(x).onRest) % Iterate over each immobility period
        tmp = [tmp; new(x).rawFP(new(x).onRest(z):new(x).offRest(z))]; % Extract and concatenate immobility periods
        waitbar(z/length(new(x).onRest),h);
    end; close(h);
    new(x).tdt_imm = tmp;
    fprintf('%s - DONE.\n',new(x).rec);
end
%% FFT
p1_mat = [];
needL = 2500*new(1).rawFs;
h = waitbar(0, 'new processing');
for x = 1:length(new)
    vec = [new(x).tdt_imm]; 
    vec = repmat(vec,[ceil(needL/length(vec)) 1]);
    vec = vec(1:needL);
    Fs = new(x).rawFs;
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
    waitbar(x/length(new),h);
end
fprintf('FFT done! \n'); close(h);

%% Normalize FFT
tmp = [];
r = [find(f == 0.1):find(f == 100)]; % Restrict to [0.01 100]
flog = log10(f(r));
f_sub = 10.^(flog); % Regenerate frequency vector from log(freq)
for x = 1:size(p1_mat,2)
    % vec_norm = normalize(log10(p1_mat(r,x)),'range'); % Normalize range from [0.01 100], scaling so range covers [0 1]
    segment = log10(p1_mat(r,x));
    vec_norm = (segment - min(segment)) / (segment(1) - min(segment));
    tmp(:,x) = vec_norm;
end
norm = tmp;

%% PLOT FFT w/o subtraction
figure; hold on
plot([flog(f == 0.5) flog(f == 0.5)],[0 1],'--k'); 
plot([flog(f == 4) flog(f == 4)],[0 1],'--k'); 
plot(flog, norm);
plot(flog, nanmean(norm_gfp,2), 'k');
xlabel('Frequency'); 
xlim([-1 1.5]); xticks([-2:2]); xticklabels({'0.01','0.1','1','10','100'});

%% PLOT FFT w/o subtraction + DOWNSAMPLE
figure; hold on
ds = 50;
plot([flog(f == 0.5) flog(f == 0.5)],[0 1],'--k'); 
plot([flog(f == 4) flog(f == 4)],[0 1],'--k'); 
shadederrbar(flog(1:ds:end), nanmean(norm((1:ds:end),:),2), SEM(norm((1:ds:end),:),2), 'r');
shadederrbar(flog(1:ds:end), nanmean(norm_gfp((1:ds:end),:),2), SEM(norm_gfp((1:ds:end),:),2), 'g');
legend({'0.5Hz','4Hz','tdTomato','','GFP',''})
xlabel('Frequency'); 
xlim([-1 1.5]); xticks([-2:2]); xticklabels({'0.01','0.1','1','10','100'});
ylabel('Power (a.u.)');
axis('square');