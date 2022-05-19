%% LOAD DATA
load('C:\Users\Anya\Desktop\IV_LOCAL\data_comb\fullcin_2022-Feb_WT.mat')
sub = cinWT; beh = behWT;
sub(1:50) = []; beh(1:15) = [];

%%
velThres = 0.25;
Fs = 50;
minRestTime = 4*Fs; minRunTime = 1*Fs;
timeThres = 4*Fs; timeShift = 0.5*Fs;

for idx_b = 1:length(beh)
    vel = beh(idx_b).vel; vel = abs(vel);
    [onsetInd,offsetInd] = getOnsetOffset(-vel,-velThres,minRunTime,minRestTime,1);
    [onsetInd,offsetInd] = adjOnsetOffset(onsetInd,offsetInd,timeThres,vel);
    onsetInd = onsetInd+timeShift; offsetInd = offsetInd-timeShift;
    beh(idx_b).onRest = onsetInd; beh(idx_b).offRest = offsetInd;
end

%% CCG rest and mvmt spikes
mat = struct; %Initialize structure to save CCG output data into
uni = unique({sub.rec}); %Find unique recording IDs across all units
Fs = 50; %Sampling frequency for behavioral data
diffFs = 1;

for x = 1:length(uni)
    fprintf('%s \n',uni{x});
    idx = find(strcmp({sub.rec},uni{x})); %Find all units that match this unique recording ID
    idx_b = find(strcmp({beh.rec},uni{x})); %Find behavior data that matches this unique recording ID
    if isempty(idx_b); continue; end
    if length(idx) < 2; continue; end %If there are <2 units (thus no possible unit pairs), continue to next unique recording ID
    ccgst = rununitccg({sub(idx).st},beh(idx_b).on/diffFs,beh(idx_b).off/diffFs,beh(idx_b).onRest/diffFs,beh(idx_b).offRest/diffFs); 

    mat(x).rec = uni{x};
    
    mat(x).n = [ccgst.pairs.n]; mat(x).m = [ccgst.pairs.m];
    mat(x).fr = [];
    mat(x).ccg = [ccgst.pairs.ccg];
    % mat(x).ccgZ = [ccgst.pairs.ccgZ];
    mat(x).shuffPrc = cell(length(ccgst.pairs),1); %5th, 50th, 95th percentile of shuffled CCG's
    
    mat(x).fr_mvmt = [];
    mat(x).ccg_mvmt = [ccgst.pairs.ccg_mvmt]; 
    % mat(x).ccgZ_mvmt = [ccgst.pairs.ccgZ_mvmt];
    mat(x).shuffPrc_mvmt = cell(length(ccgst.pairs),1); %5th, 50th, 95th percentile of shuffled CCG's
    
    mat(x).fr_rest = [];
    mat(x).ccg_rest = [ccgst.pairs.ccg_rest]; 
    % mat(x).ccgZ_rest = [ccgst.pairs.ccgZ_rest];
    mat(x).shuffPrc_rest = cell(length(ccgst.pairs),1); %5th, 50th, 95th percentile of shuffled CCG's
    
    mat(x).dist = [];
    for y = 1:length(ccgst.pairs)
        mat(x).fr(y) = 1/mean(diff(ccgst.times(ccgst.pairs(y).m).full));
        mat(x).fr_mvmt(y) = 1/mean(diff(extractEventST(ccgst.times(ccgst.pairs(y).m).full,beh(idx_b).on/diffFs,beh(idx_b).off/diffFs,0)));
        mat(x).fr_rest(y) = 1/mean(diff(extractEventST(ccgst.times(ccgst.pairs(y).m).full,beh(idx_b).onRest/diffFs,beh(idx_b).offRest/diffFs,0)));
        mat(x).shuffPrc{y} = prctile(ccgst.pairs(y).shuff,[2.5 50 97.5],2); %5th, 50th, 95th percentile of shuffled CCG's
        mat(x).shuffPrc_mvmt{y} = prctile(ccgst.pairs(y).shuff_mvmt,[2.5 50 97.5],2); %5th, 50th, 95th percentile of shuffled CCG's
        mat(x).shuffPrc_rest{y} = prctile(ccgst.pairs(y).shuff_rest,[2.5 50 97.5],2); %5th, 50th, 95th percentile of shuffled CCG's
%         xc_diff = sub(idx(ccgst.pairs(y).n)).coor(1) - sub(idx(ccgst.pairs(y).m)).coor(1); %Distance in x- or y-dimension
%         zc_diff = sub(idx(ccgst.pairs(y).n)).coor(2) - sub(idx(ccgst.pairs(y).m)).coor(2); %Distnace in z-dimension
%         mat(x).dist(y) = hypot(xc_diff,zc_diff);
    end
end
time = ccgst.lag;
