%%
% load('C:\Users\Anya\Desktop\FP_LOCAL\ach+da\FFT_pulldata_fName+fPath.mat'); % CONTROL DATA
load('C:\Users\Anya\Desktop\FP_LOCAL\cannula\raw_data\FFT_ACh+scop+ChAT+GFP_normForPlot.mat');
load('C:\Users\Anya\Desktop\FP_LOCAL\lesion\da_rawfp_DA_IMM_normFFToutput.mat');
%%
fPath = 'C:\Users\Anya\Desktop\FP_LOCAL\'; 
[fName,fPath] = uigetfile([fPath,'*.mat'],'MultiSelect','On');
if ~iscell(fName); fName = {fName}; end

%%
new = struct;
h = waitbar(0, 'new processing');
for y = 1:length(fName)
    load(fullfile(fPath,fName{y})); % Load data file
    [an,b] = strtok(fName{y},'_'); day = strtok(b,'_'); % Parse file name
    x = 1 + length(new);
    new(x).rec = [an,'-',day]; 
    new(x).site = 'DLS';
    new(x).rx = 'control';
    
    %% Pull parameters required for this analysis
    if isfield(data.gen,'params')
        params = data.gen.params; % Extract params structure
        dsRate = params.dsRate; dsType = params.dsType; % General downsampling parameter
        rawFs = data.gen.acqFs; Fs = data.gen.Fs;
    end
    %% Process photometry data
    fprintf('Processing %s ... ',new(x).rec);
    new(x).FPnames = data.acq.FPnames;
    new(x).rawFP = data.acq.FP;
    fprintf('DONE.\n');
    new(x).rawFs = rawFs;
    if isfield(data,'final'); if isfield(data.final,'mov')
    new(x).on = data.final.mov.onsets*dsRate;
    new(x).off = data.final.mov.offsets*dsRate;
        new(x).onRest = data.final.mov.onsetsRest*dsRate;
        new(x).offRest = data.final.mov.offsetsRest*dsRate;
        end; end
%     if isfield(data.final,'rew'); if isfield(data.final.rew,'delivery')
%         new(x).reward = data.final.rew.delivery*dsRate;
%         end; end
    waitbar(y/length(fName),h);
end
close(h);
if isempty(new(1).rawFs); new(1) = []; end

%% extract during IMMOBILITY + during INFUSION window
for x = 1:length(new)
    tmp = [];
    for z = 1:length(new(x).onRest) % Iterate over each immobility period
        tmp = [tmp; new(x).rawFP{1}(new(x).onRest(z):new(x).offRest(z))]; % Extract and concatenate immobility periods
%         for z = 1:length(new(x).on)
%         tmp = [tmp; new(x).rawFP{1}(new(x).on(z):new(x).off(z))]; % Extract and concatenate immobility periods
    end
    new(x).ach_imm = tmp;
end
fprintf('Extraction done! \n');

%% FFT
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
    tmp(:,x) = vec_norm;
end
norm_da = tmp;
fprintf('Normalization done! \n');

%% Subtract gfp
sub_1 = [norm_da_imm, norm_da_mov]; % SUBTRACT
sub_2 = norm_d2r; % GFP or tdTomato

sub_gfp = []; for x = 1:size(sub_1,2); sub_gfp(:,x) = sub_1(:,x) - nanmean(sub_2,2); end % Subtract avg FFT for mAChR antagonist

%% AUC post-subtraction
auc = [];
f_sub = 10.^(flog); % Regenerate frequency vector from log(freq)
r_14 = [find(f_sub == 0.5):find(f_sub == 4)]; % AUC from [1 4] Hz
for x = 1:size(sub_gfp,2)
    auc(x) = trapz(sub_gfp(r_14,x))/length(r_14);
end

%% STATISTICS - FIG 1
sub_1 = [norm_ach_an(:,[6:15])]; % SUBTRACT
sub_2 = norm_gfp_an; % GFP or tdTomato
sub_gfp = []; for x = 1:size(sub_1,2); sub_gfp(:,x) = sub_1(:,x) - nanmean(sub_2,2); end % Subtract avg FFT for mAChR antagonist
[a,b] = max(sub_gfp);
fprintf('mode ACh: %1.2f +/- %1.2f Hz \n',nanmean(f(b)),SEM(f(b),2))

sub_1 = [norm_da]; % SUBTRACT
sub_2 = norm_tdt; % GFP or tdTomato
sub_gfp = []; for x = 1:size(sub_1,2); sub_gfp(:,x) = sub_1(:,x) - nanmean(sub_2,2); end % Subtract avg FFT for mAChR antagonist
[a,b] = max(sub_gfp);
fprintf('mode DA: %1.2f +/- %1.2f Hz \n',nanmean(f(b)),SEM(f(b),2))

sub_1 = [norm_ach_an(:,[6:15]), norm_scop, norm_chat_an]; % SUBTRACT
sub_2 = norm_gfp_an; % GFP or tdTomato
sub_gfp = []; for x = 1:size(sub_1,2); sub_gfp(:,x) = sub_1(:,x) - nanmean(sub_2,2); end % Subtract avg FFT for mAChR antagonis
auc = [];
f_sub = 10.^(flog); % Regenerate frequency vector from log(freq)
r_14 = [find(f_sub == 0.5):find(f_sub == 4)]; % AUC from [1 4] Hz
for x = 1:size(sub_gfp,2)
    auc(x) = trapz(sub_gfp(r_14,x))/length(r_14);
end
auc(auc < 0) = 0;
[p,~,stats] = kruskalwallis(auc', [ones(10,1);2*ones(5,1);3*ones(5,1)], 'off');
c = multcompare(stats,'display','off');
fprintf('ACh vs antag: %1.4f, ACh vs ChATcKO: %1.4f \n',c(1,6),c(2,6))
fprintf('ranksum ACh vs antag: %1.4f, ACh vs ChATcKO: %1.4f \n',ranksum(auc(1:10),auc(11:15)),ranksum(auc(1:10),auc(16:20)))
[~,p2] = ttest2(auc(1:10),auc(11:15)); [~,p2(2)] = ttest2(auc(1:10),auc(16:20));
fprintf('ttest2 ACh vs antag: %1.4f, ACh vs ChATcKO: %1.4f \n',p2(1),p2(2))

sub_1 = [norm_da, norm_d2r, norm_lesion]; % SUBTRACT
sub_2 = norm_d2r; % GFP or tdTomato
sub_gfp = []; for x = 1:size(sub_1,2); sub_gfp(:,x) = sub_1(:,x) - nanmean(sub_2,2); end % Subtract avg FFT for mAChR antagonis
auc = [];
f_sub = 10.^(flog); % Regenerate frequency vector from log(freq)
r_14 = [find(f_sub == 0.5):find(f_sub == 4)]; % AUC from [1 4] Hz
for x = 1:size(sub_gfp,2)
    auc(x) = trapz(sub_gfp(r_14,x))/length(r_14);
end
auc(auc < 0) = 0;
[p,~,stats] = kruskalwallis(auc', [ones(10,1);2*ones(8,1);3*ones(6,1)], 'off');
c = multcompare(stats,'display','off');
fprintf('DA vs antag: %1.4f, DA vs lesion: %1.4f \n',c(1,6),c(2,6))
fprintf('ranksum DA vs antag: %1.4f, DA vs lesion: %1.4f \n',ranksum(auc(1:10),auc(11:18)),ranksum(auc(1:10),auc(19:24)))
[~,p2] = ttest2(auc(1:10),auc(11:18)); [~,p2(2)] = ttest2(auc(1:10),auc(19:24));
fprintf('ttest2 DA vs antag: %1.4f, DA vs lesion: %1.4f \n',p2(1),p2(2))