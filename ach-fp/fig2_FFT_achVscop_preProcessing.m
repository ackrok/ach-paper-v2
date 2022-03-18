% subtract FFT
% Processing: digitalLIA, but commented out all filtering (lines 74, 78, 80)
% Anya Krok, Jan 25 2022

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
    new(x).site = 'DMS';
    new(x).rx = '';
    
    %% Pull parameters required for this analysis
params = data.gen.params; % Extract params structure
%     lpCut = 500; filtOrder = params.FP.filtOrder; % Filter Properties
dsRate = params.dsRate; dsType = params.dsType; % General downsampling parameter
%     sigEdge = params.FP.sigEdge; % Demodulation-specific property
rawFs = data.gen.acqFs; Fs = rawFs/dsRate;
%     nFP = length(data.acq.FP); %Obtain number of FP channels (Includes red control)
    refSig = data.acq.refSig; %Pull reference signals
    
    %% Process photometry data
    fprintf('Processing %s ... ',new(x).rec);
    new(x).rawFP = data.acq.FP;
%     demod_dual = cell(nFP,1);
%     for z = 1
%         fprintf('nFP %d of %d: ',z,nFP);
%         rawFP = data.acq.FP{z}; %Extract FP trace
%         modFreq = params.FP.modFreq(z); % ANYA EDIT 21/04/06
%         ref = findRef(modFreq,refSig,rawFs); %Find the reference signal from the refsig array using modulation frequency
%         fprintf('ref sig found ... ');
%         % NOTE: digitalLIA comment out all filtering (lines 74, 78, 80)
%         demod = digitalLIA(rawFP,ref,modFreq,rawFs,lpCut,filtOrder); %Peform the demodulation
%         if sigEdge ~= 0 %Remove the beginning and the edge if the setting isn't 0
%             demod = demod((sigEdge*rawFs)+1:end-(sigEdge*rawFs));
%         end
%         fprintf('demod complete ... ');
%         demod_dual{z} = demod;
%     end
    fprintf('DONE.\n');
%     new(x).demod = demod_dual{z};
    new(x).rawFs = rawFs;
%     new(x).on = data.final.mov.onsets*dsRate;
%     new(x).off = data.final.mov.offsets*dsRate;
    if isfield(data.final,'mov')
        new(x).onRest = data.final.mov.onsetsRest*dsRate;
        new(x).offRest = data.final.mov.offsetsRest*dsRate;
    end
%     if isfield(data.final,'rew'); if isfield(data.final.rew,'delivery')
%         new(x).reward = data.final.rew.delivery*dsRate;
%         end; end
    waitbar(y/length(fName),h);
end
close(h);
if isempty(new(1).rawFs); new(1) = []; end

%%
for x = [3 5 9 11]; new(x).rx = 'aCSF'; end
for x = [2 6 8 12]; new(x).rx = 'd1d2'; end
for x = [1 4 7 10]; new(x).rx = 'glu'; end
new = new([3 5 9 11 2 6 8 12 1 4 7 10]);

%% separate ACh, DA data to make structure smaller
new_ach = struct; new_da = struct;
for x = 1:length(new)
    new_ach(x).rec = new(x).rec;
    new_ach(x).site = new(x).site;
    new_ach(x).rx = new(x).rx;
    new_ach(x).demod = new(x).demod{1};
    new_ach(x).rawFs = new(x).rawFs;
    
    new_da(x).rec = new(x).rec;
    new_da(x).site = new(x).site;
    new_da(x).rx = new(x).rx;
    new_da(x).demod = new(x).demod{2};
    new_da(x).rawFs = new(x).rawFs;
end