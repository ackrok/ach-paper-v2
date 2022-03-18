%%
for x = 1:length(lfp_can)
    lfp_imm = [];
    Fs = lfp_can(x).Fs;
    for z = 1:length(lfp_can(x).onRest)
        lfp_imm = [lfp_imm; lfp_can(x).lfp(lfp_can(x).onRest(z)*(Fs/50):lfp_can(x).offRest(z)*(Fs/50))]; % Extract and concatenate immobility periods
    end
    lfp_can(x).lfp_imm = lfp_imm;
end
fprintf('Done \n');

%%
p1_mat = [];
for x = 1:length(lfp_can)
    vec = lfp_can(x).lfp_imm;
    Fs = lfp_can(x).Fs;

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
    p1_mat = [p1_mat, movmean(P1,500)];
end
fprintf('FFT done \n');

%% PLOT RAW
figure; hold on
ds = 50;
% plot([flog(f == 0.5) flog(f == 0.5)],[0 1],'--k'); 
% plot([flog(f == 4) flog(f == 4)],[0 1],'--k'); 
plot(f(1:ds:end), p1_mat(1:ds:end,1), 'g');
plot(f(1:ds:end), p1_mat(1:ds:end,2), 'r');
plot(f(1:ds:end), p1_mat(1:ds:end,3), 'c');
plot(f(1:ds:end), p1_mat(1:ds:end,4), 'm');
plot(f(1:ds:end), p1_mat(1:ds:end,5), 'b');
plot([0.5 0.5],[0 2],'--k'); plot([4 4],[0 2],'--k');
legend({'aCSF-1','glu-1','aCSF-2','glu-2','aCSF2b','0.5Hz','4Hz'})
xlabel('Frequency'); xlim([0 25]);
ylabel('Power (a.u.)');
title('cannula LFP'); axis('square');

%%
tmp = [];
r = [find(f == 0.1):find(f == 100)]; % Restrict to [0.01 100]
flog = log10(f(r));
f_sub = 10.^(flog); % Regenerate frequency vector from log(freq)
for x = 1:size(p1_mat,2)
    vec_norm = normalize(log10(p1_mat(r,x)),'range'); % Normalize range from [0.01 100], scaling so range covers [0 1]
    tmp(:,x) = vec_norm;
end
norm = tmp;
fprintf('Normalization done \n');

%% Plot NORMALIZATION
figure; hold on
ds = 1;
plot([flog(f == 0.5) flog(f == 0.5)],[0 1],'--k'); 
plot([flog(f == 4) flog(f == 4)],[0 1],'--k'); 
plot(flog(1:ds:end), norm(1:ds:end,1), 'g');
plot(flog(1:ds:end), norm(1:ds:end,2), 'r');
plot(flog(1:ds:end), norm(1:ds:end,3), 'c');
plot(flog(1:ds:end), norm(1:ds:end,4), 'm');
plot(flog(1:ds:end), norm(1:ds:end,5), 'b');
legend({'0.5Hz','4Hz','aCSF','glu','aCSF-2','glu-2'})
xlabel('Frequency'); 
xlim([-1 flog(f == 100)]); xticks([-2:2]); xticklabels({'0.01','0.1','1','10','100'});
ylabel('Power (a.u.)');
title('cannula LFP');
axis('square');