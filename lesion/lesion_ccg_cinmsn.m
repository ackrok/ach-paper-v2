%% CCG - CIN vs MSN
sub = cinLes;
msn = lesion(find([lesion.label] == 1));

%%
mat = struct; %Initialize structure to save CCG output data into
uni = unique({sub.rec}); %Find unique recording IDs across all units
Fs = 50; %Sampling frequency for behavioral data
diffFs = 50;

for x = 1:length(uni)
    fprintf('%s \n',uni{x});
    idx = find(strcmp({sub.rec},uni{x})); %Find all units that match this unique recording ID
    idx_msn = find(strcmp({msn.rec},uni{x})); 
    if isempty(idx_msn); continue; end
    idx_b = find(strcmp({beh.rec},uni{x})); %Find behavior data that matches this unique recording ID
    if isempty(idx_b); continue; end
    if length(idx) < 2; continue; end %If there are <2 units (thus no possible unit pairs), continue to next unique recording ID
    ccgst = rununitccg_cin({msn(idx_msn).st},{sub(idx).st},beh(idx_b).on/diffFs,beh(idx_b).off/diffFs,beh(idx_b).onRest/diffFs,beh(idx_b).offRest/diffFs); 
%%
    for y = 1:length(ccgst.pairs)
        n = y;%1+length(mat);
        mat(n).rec = uni{x};
        mat(n).n = ccgst.pairs(y).n;
        mat(n).m = ccgst.pairs(y).m;
        mat(n).fr = []; mat(n).fr_mvmt = []; mat(n).fr_rest = [];
        mat(n).ccg = ccgst.pairs(y).ccg;
        mat(n).ccg_mvmt = ccgst.pairs(y).ccg_mvmt;
        mat(n).ccg_rest = ccgst.pairs(y).ccg_rest;
        for z = 1:length(idx_msn)
            mat(n).fr(z) = 1/mean(diff(ccgst.times(z).full));
            mat(n).fr_mvmt(z) = 1/mean(diff(extractEventST(ccgst.times(z).full,beh(idx_b).on/diffFs,beh(idx_b).off/diffFs,0)));
            mat(n).fr_rest(z) = 1/mean(diff(extractEventST(ccgst.times(z).full,beh(idx_b).onRest/diffFs,beh(idx_b).offRest/diffFs,0)));
        end
        mat(n).ccgDelta = (mat(n).ccg - mat(n).fr)./mat(n).fr;
        mat(n).ccgDelta_mvmt = (mat(n).ccg_mvmt - mat(n).fr_mvmt)./mat(n).fr_mvmt;
        mat(n).ccgDelta_rest = (mat(n).ccg_rest - mat(n).fr_rest)./mat(n).fr_rest;
        mat(n).shuffPrc = prctile(ccgst.pairs(y).shuff,[5 50 95],2); %5th, 50th, 95th percentile of shuffled CCG's
        mat(n).shuffDelta = (mat(n).shuffPrc - nanmean(mat(n).fr))./nanmean(mat(n).fr);
    end
end
time = ccgst.lag;

%%
x = 1;
sm = 10;

figure;
shadederrbar(time, movmean(nanmean(mat(x).ccgDelta_rest,2),sm), movmean(SEM(mat(x).ccgDelta_rest,2),sm), 'b');
