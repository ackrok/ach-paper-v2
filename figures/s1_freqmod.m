sampling_rate = 5000;
time = 0:1/sampling_rate:0.02;
freq = 217; s1 =  sin(2*pi*time*freq);
freq = 319; s2 =  sin(2*pi*time*freq);
%general formula : Amplitude*sin(2*pi*freq*time)
figure(1),clf; hold on;
plot(time.*1000,s1,'g'); plot(time.*1000,s2,'m');
xlabel('time (ms)')
title('Sine Wave')

%%
x = 39;
load('C:\Users\Anya\Desktop\FP_LOCAL\ach+da\AK190_210810_0001.mat'); % load raw data
%%
win = [7 8]; % window for plotting
r = [find(data.acq.time == win(1)):find(data.acq.time == win(2))]; % segment for raw data
r2 = [find(data.final.time == win(1)):find(data.final.time == win(2))]; % segment for processed data
t = [0:1/5000:diff(win)]; t2 = [0:1/50:diff(win)];

win = [7 7.05];
r3 = [find(data.acq.time == win(1)):find(data.acq.time == win(2))]; % segment for raw data
t3 = [0:1/5000:0.05]; 
clr = {'g','m'};

fig = figure;
subplot(3,2,1); hold on
for y = 1:2; plot(t, normalize(data.acq.refSig{y}(r),'range'), clr{y}); end
title('reference signal');

subplot(3,2,2); hold on
for y = 1:2; plot(t3, normalize(data.acq.refSig{y}(r3),'range'), clr{y}); end
title('zoom in: reference signal');

subplot(3,2,3); hold on
for y = 1:2; plot(t, data.acq.FP{y}(r), clr{y}); end
title('raw data');

subplot(3,2,4); hold on
for y = 1:2; plot(t3, data.acq.FP{y}(r3), clr{y}); end
title('zoom in: raw data');

subplot(3,2,5); hold on
for y = 1:2; plot(t2, data.final.FP{y}(r2), clr{y}); end
title('raw data');