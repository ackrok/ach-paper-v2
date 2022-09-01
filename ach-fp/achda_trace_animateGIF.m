x = 36; %JM007 210610
 
fp_1 = beh(x).FP{1}; Fs = beh(x).Fs; % extract photometry signal from structure
fp_mu = nanmean(fp_1); % mean of entire photometry signal
fp_1 = fp_1 - fp_mu; % SUBTRACT baseline% from fp%, now centered on BASELINE
fp_2 = beh(x).FP{2}; % dopamine
fp_2 = fp_2 - nanmean(fp_2);
f = [0.5 4];
fp_1_filt = filterFP(fp_1,Fs,f,10,'bandpass');
fp_2_filt = filterFP(fp_2,Fs,f,10,'bandpass');

win = [277 279]; % x = 30, immobility
% win = [505 507]; % x = 30, locomotion
% win = [1441.5 1442.6]; % x = 30, locomotion
% win = [785 786]; % x = 15, reward
win = win.*Fs; nSamp = win(2) - win(1) + 1;

win_fp = [];
win_fp(:,1) = fp_1_filt(win(1):win(2));
win_fp(:,2) = fp_2_filt(win(1):win(2));

win_fp_norm = [];
win_fp_norm(:,1) = normalize(win_fp(:,1),'range');
win_fp_norm(:,2) = normalize(win_fp(:,2),'range');
win_fp_norm = win_fp_norm - nanmean(win_fp_norm);

%% TEST
figure;
clr = jet(nSamp);
subplot(2,1,2); 
plot(win_fp_norm(:,1),'k', 'Color', [0 0 0]); hold on
xlim([0 100]); xticks([0 50 100]); xticklabels({'0','1','2'});
ylabel('ACh'); yticks([-1 1]);
subplot(2,1,1); 
plot(win_fp_norm(:,2),'k', 'Color', [0 0 0]); hold on
xlim([0 100]); xticks([0 50 100]); xticklabels({'0','1','2'});
ylabel('DA'); yticks([-1 1]);
set(gcf,'color','w');
for z = 1:nSamp
    subplot(2,1,1); 
    plot(z,win_fp_norm(z,1),'.','MarkerSize',20,'Color',clr(z,:));
    subplot(2,1,2); 
    plot(z,win_fp_norm(z,2),'.','MarkerSize',20,'Color',clr(z,:));
end

%% both plotted
figure(1)

filename = 'C:\Users\Anya\Desktop\FP_LOCAL\achda_trace_GIF.gif';
clr = jet(nSamp);
% subplot(2,1,1); 
    plot(win_fp_norm(:,2),'k', 'Color', 'k', 'LineWidth', 1.5); hold on
    xlim([0 100]); xticks([0 50 100]); xticklabels({'0','1','2'});
    ylabel('DA'); yticks([-1 1]);
    title('ACh');
% subplot(2,1,2); 
%     plot(win_fp_norm(:,1),'k', 'Color', [0 0 0]); hold on
%     xlim([0 100]); xticks([0 50 100]); xticklabels({'0','1','2'});
%     ylabel('ACh'); yticks([-1 1]);
set(gcf,'color','w');

for z = 1:nSamp
%     subplot(2,1,1); 
        plot(z,win_fp_norm(z,1),'.','MarkerSize',30,'Color',clr(z,:));
%     subplot(2,1,2); 
%         plot(z,win_fp_norm(z,2),'.','MarkerSize',20,'Color',clr(z,:));
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if z == 1
        imwrite(imind,cm,filename,'gif', 'DelayTime', 1/30, 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif', 'DelayTime', 1/30, 'WriteMode','append');
    end
end

%% both plotted
figure(1)

filename = 'C:\Users\Anya\Desktop\FP_LOCAL\achda_trace_GIF.gif';
clr = jet(nSamp);
subplot(2,1,1); 
    plot(win_fp_norm(:,2),'k', 'Color', [0 0 0]); hold on
    xlim([0 100]); xticks([0 50 100]); xticklabels({'0','1','2'});
    ylabel('DA'); yticks([-1 1]);
subplot(2,1,2); 
    plot(win_fp_norm(:,1),'k', 'Color', [0 0 0]); hold on
    xlim([0 100]); xticks([0 50 100]); xticklabels({'0','1','2'});
    ylabel('ACh'); yticks([-1 1]);
set(gcf,'color','w');

for z = 1:nSamp
    subplot(2,1,1); 
        plot(z,win_fp_norm(z,1),'.','MarkerSize',20,'Color',clr(z,:));
    subplot(2,1,2); 
        plot(z,win_fp_norm(z,2),'.','MarkerSize',20,'Color',clr(z,:));
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if z == 1
        imwrite(imind,cm,filename,'gif', 'DelayTime', 1/30, 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif', 'DelayTime', 1/30, 'WriteMode','append');
    end
end

%% phase plot
signal = fp_1;
Fpass = [0.5 4];
Fs = 50; %sampling rate, has to be at least double of your high pass frequency
Wn = [Fpass(1)/(Fs/2) Fpass(2)/(Fs/2)]; 
[b,a] = butter(3,Wn);
data_filt= filtfilt(b,a,signal); % signal is your lfp data, output is the filtered data
H = hilbert(double(data_filt));
data_phase = angle(H); % output is the instantaneous phase
fp_phase = data_phase;
fp_deg = rad2deg(data_phase);
win_fp_deg = fp_deg(win(1):win(2));

figure;
plot(win_fp_norm(:,2),'k', 'Color', [0 0 0]); hold on
plot(win_fp_norm(:,1),'k', 'Color', [0 0 0 0.2]); hold on
clr = jet(length([-180:180]));
for z = 1:nSamp
    plot(z,win_fp_norm(z,1),'.','MarkerSize',20,'Color',clr(181+round(win_fp_deg(z)),:));
end
xlim([0 100]); xticks([0 50 100]); xticklabels({'0','1','2'});
ylabel('DA'); yticks([-1 1]);
set(gcf,'color','w');