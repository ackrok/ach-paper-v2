[fName,fPath] = uigetfile('*.mat','MultiSelect','On');
%%
norm = [];
beh = behwt; behFs = 50;

for x = 1:length(fName)
    load(fullfile(fPath,fName{x}));
    
    lfp_mat = double(lfp.data);
    Fs = lfp.samplingRate;

    lfp_sub = [];
    b = find(strcmp({beh.rec},lfp.rec));
    imm = [beh(b).onRest(:),beh(b).offRest(:)].*(Fs/behFs);
    for z = 1:size(imm,1)
        lfp_sub = [lfp_sub; lfp_mat(imm(z,1):imm(z,2),:)];
    end

    %% FFT
    p1_mat = [];

    for y = 1:size(lfp_mat,2)
    vec = lfp_sub(:,y); 
    Fs = lfp.samplingRate;

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
    p1_mat(:,y) = [movmean(P1,500)];

    end
    fprintf('FFT done! \n');

%     figure; plot(f, p1_mat);
%     xlim([0 100]);
%     ylabel('power'); xlabel('frequency');

    %% Normalize FFT
    tmp = [];
    r = [find(f == 0.1):find(f == 100)]; % Restrict to [0.01 100]
    flog = log10(f(r));
    f_sub = 10.^(flog); % Regenerate frequency vector from log(freq)
    for x = 1:size(p1_mat,2)
        a = log10(p1_mat(r,x));
        vec_norm = (a - a(end))./(a(1) - a(end)); 
        tmp(:,x) = vec_norm;
    end
    norm = [norm, tmp];
    fprintf('Normalization done! \n');

%     figure; plot(flog, norm);
%     %xlim([0 100]);
%     xlabel('Frequency'); 
%     xlim([-1 2]); xticks([-2:2]); xticklabels({'0.01','0.1','1','10','100'});
%     ylabel('Power (a.u.)');
    
end

%%
figure; hold on
plot(flog, nanmean(norm,2))
%shadederrbar(flog, nanmean(norm,2), SEM(norm,2), 'r');
plot([flog(f == 0.5) flog(f == 0.5)],[0 5],'--k'); 
plot([flog(f == 4) flog(f == 4)],[0 5],'--k'); 
xlabel('Frequency'); 
xlim([-1 2]); xticks([-2:2]); xticklabels({'0.01','0.1','1','10','100'});
ylabel('Power (a.u.)');