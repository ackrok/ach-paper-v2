% I actually think the attached function should do the trick for you.
% You give it spikes (cell array, each entry holds spike times for a given neuron), and lfp... Here it's expecting formatting
% that's provided by the bz_getLFP function in buzcode. Let me know if you run into issues with this.
% You give a passband for your frequencies of interest (mine are [6 14] for hippocampal theta), and you give a powerthreshold.
% This is a threshold on the zscored power in the frequency band - you will only be considering phase modulation of your spikes
% when there is sufficient power. The default is a threshold of 2, which works well.
% 
% Some things about the output:
% phasestats.m - holds neuron's preferred phase (in radians)
% phasestats.p - the p-value for the phase modulation (Rayleigh test)
% phasedistros - percentage distribution of spikes binned across phases (this will be useful for your example)
% phasebins - the phase bins, your x axis ; it's in radians, 0 to 2pi, with 0 being the peak ; if you want to plot a sinusoid to illustrate the LFP, plot the cosine function (i.e., cos( phasebins ) )
% 
% We usually show two cycles when displaying an example neuron's phase locking (I'm making the x axis into degrees, rather than radians):
% plot([rad2deg(phasebins ) 360+rad2deg (phasebins )], [ phasedistros(:,n)    phasedistros(:,n)  ]), where n is the neuron of interest (note - this line is just to give you the idea, the syntax might be a little off) 
% hold on
% plot([rad2deg(phasebins ) 360+rad2deg (phasebins )],  [ cos(phasebins )  cos(phasebins )]  % show the cosine function to illustrate the LFP (you can also find example lfp for this, but that is a tedious endeavor)
% 
% If you want to show the phase distribution across all cells, collect phasestats.m across all cells, and do a histogram using phasebins as the bins ; again, if you want everything to be in degrees, 
% you can just convert using rad2deg as is illustrated above

rec = 'IV055_rec01';

lfp = struct; ch = [39];
lfp.data = bz_LoadBinary(['C:\Users\Anya\Desktop\IV_LOCAL\ACh\',rec,'.dat'],...
        'duration',Inf,'start',0,'frequency',30000,...
        'nchannels',128,'channels',ch,...
                  'downsample',24);
lfp.channels = ch;
lfp.samplingRate = 1250;
lfp.rec = rec;
fprintf('LFP loaded.\n');

save(['C:\Users\Anya\Desktop\IV_LOCAL\ACh\',rec,'_LFP_ch',num2str(ch),'.mat'],'lfp');

%%
band = [0.5 4];
order = 10; %10 or 12 
filt = designfilt('bandpassiir','FilterOrder',order,'HalfPowerFrequency1',band(1),'HalfPowerFrequency2',band(2),'SampleRate',lfp.samplingRate,'DesignMethod','butter');
lfp_band = [];
for x = 1:2
    lfp_band(:,x) = filtfilt(filt,double(lfp.data(:,x)));
end

figure; hold on
plot(lfp_band(:,1)); plot(lfp_band(:,2));
[corr,lags] = xcorr(lfp_band(:,1), lfp_band(:,2), 1250, 'coeff');
figure; plot(lags/1250, corr);


%%
st = {cinWT(strcmp({cinWT.rec},lfp.rec)).st};
spikes = struct; spikes.times = st;

sub = WT(find([WT.label] == 1));
spikes = struct; spikes.times = {sub(strcmp({sub.rec},lfp.rec)).st};

%%
data = bz_PhaseModulation('spikes',spikes, 'lfp',lfp, 'passband', [0.5 4], ...
    'intervals',[0 inf], 'samplingRate',1250, 'method','hilbert', ...
    'numBins',30, 'powerThresh',1, 'plotting',false, 'saveMat',false);
fprintf('Spikes - LFP done.\n');

%%
fp = struct;
fp.data = behACh(strcmp({behACh.rec},rec)).FP{1};
fp.channels = 1;
fp.samplingRate = 50;
st = spikes.times;
test = cell(length(st),1);
for x = 1:length(st)
    idx = find(st{x}.*fp.samplingRate > length(fp.data),1);
    if ~isempty(idx); test{x} = st{x}(1:idx-1); else test{x} = st{x}; end
end
spikes.times = test;

data = bz_PhaseModulation('spikes',spikes, 'lfp',fp, 'passband', [0.1 10], ...
    'intervals',[0 inf], 'samplingRate',50, 'method','hilbert', ...
    'numBins',30, 'powerThresh',0.5, 'plotting',false, 'saveMat',false);
fprintf('Spikes - FP done.\n');

%%
figure;
phasebins = data.phasebins;
phasedistros = data.phasedistros;
for n = 1:length(spikes.times)
    sp(n) = subplot(6,6,n);
    plot([rad2deg(phasebins) 360+rad2deg(phasebins)], [phasedistros(:,n) phasedistros(:,n)]);
    yyaxis right
    plot([rad2deg(phasebins) 360+rad2deg(phasebins)],  [cos(phasebins) cos(phasebins)]);    
    title(sprintf('%d: phase = %d, p = %1.2f',n,round(rad2deg(data.phasestats.m(n))),data.phasestats.p(n)));
    xlim([0 720]); xticks([0:180:720]);
end