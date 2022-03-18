[fName_raw, fPath_raw] = uigetfile('*.mat','Select RAW data files','MultiSelect','On');
[fName_fin, fPath_fin] = uigetfile('*.mat','Select processed data files','MultiSelect','On');

%%
params_fp = struct;
params_fp.lpCut = 50; params_fp.filtOrder = 8;
% params.rawFs = 2000; params.dsRate = 40; params.dsType = 2;
params_fp.interpType = 'linear'; params_fp.fitType = 'interp'; 
params_fp.basePrc = 5; params_fp.winSize = 10; params_fp.winOv = 0;

for x = 1:length(fName_raw)
    fprintf('%s: ',fName_raw{x});
    load(fullfile(fPath_raw,fName_raw{x})); fprintf('dataRaw loaded. ');
    idx = find(strcmp(fName_fin,[dataRaw.mouserec,'.mat']));
    load(fullfile(fPath_fin,fName_fin{idx})); 
    
    data.gen.params.FP = params_fp;
    data.gen.params.dsType = 2;
    [dataRaw] = processFP(dataRaw, data.gen.params); fprintf('processFP done. ');
    data.final.FP = dataRaw.final.FP;
    data.final.nbFP = dataRaw.final.nbFP;
    data.final.FPbaseline = dataRaw.final.FPbaseline;
    
    data.gen.params.mov.velThres_rest = 0.25;
    data.gen.params.mov.minRestTime_rest = 4; 
    data.gen.params.mov.minRunTime_rest = 1;
    data.gen.params.mov.timeThres_rest = 4; 
    data.gen.params.mov.timeShift_rest = 0.5;
    data = processOnsetOffset(data, data.gen.params);
    data = processRestOnsetOffset(data, data.gen.params);
    
    save(fullfile(fPath_fin,fName_fin{idx}),'data'); fprintf('saved.\n');
end

%%
behDA = struct;
for x = 1:length(fName_fin)
    load(fullfile(fPath_fin,fName_fin{x}));
    y = 1 + length(behDA);
    behDA(y).rec = data.mouserec;
    behDA(y).Fs = data.gen.Fs;
    behDA(y).time = data.final.time;
    behDA(y).FP = data.final.FP;
    behDA(y).nbFP = data.final.nbFP;
    behDA(y).FPnames = data.final.FPnames;
    behDA(y).vel = data.final.vel;
    behDA(y).on = data.final.mov.onsets;
    behDA(y).off = data.final.mov.offsets;
    if data.gen.params.mov.velThres_rest ~= 0.25
        data.gen.params.mov.velThres_rest = 0.25;
        data.gen.params.mov.minRestTime_rest = 4; 
        data.gen.params.mov.minRunTime_rest = 1;
        data.gen.params.mov.timeThres_rest = 4; 
        data.gen.params.mov.timeShift_rest = 0.5;
        data = processRestOnsetOffset(data, data.gen.params);
        save(fullfile(fPath_fin,fName_fin{x}),'data'); 
    end
    behDA(y).onRest = data.final.mov.onsetsRest;
    behDA(y).offRest = data.final.mov.offsetsRest;
end
if isempty(behDA(1).Fs); behDA(1) = []; end