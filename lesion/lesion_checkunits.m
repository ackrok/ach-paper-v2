figure;
for x = 1:length(cinLes)
    sp(x) = subplot(7,7,x);
    plot(cinLes(x).wf);
    title(sprintf('%s %d',cinLes(x).rec,cinLes(x).n));
end

%%
figure; violinplot([cinLes.fr]);

%%
figure; plot(nanmean(ccgDelta_rest,2));

%%
figure;
for x = 1:length(mat)
    sp(x) = subplot(5,3,x); hold on
    tmp = (mat(x).ccg_rest - mat(x).fr_rest)./mat(x).fr_rest;
    plot(time, tmp); xlim([-0.5 0.5]);
    title(sprintf('%s',mat(x).rec));
end

%%
x = 18;
figure;
unitST = cinLes(x).st; nST = length(unitST); %Extract spike times for this unit
binw = 20; %Default: bin width = 5 seconds
[~, edges] = histcounts(unitST,'BinWidth',binw); %Edges of histogram based on specified bin width 
h1 = histogram(unitST, edges, 'Normalization','countdensity'); %Histogram of spike times binned in Xs bins, as specified by binw variable
h1.FaceColor = [0.95 0.95 0.95]; 

%%
figure;
for x = 1:6
    subplot(3,2,x);
    % plot(cinLes(x+18).wf);
    unitST = cinLes(x+18).st; nST = length(unitST); %Extract spike times for this unit
    binw = 20; %Default: bin width = 5 seconds
    [~, edges] = histcounts(unitST,'BinWidth',binw); %Edges of histogram based on specified bin width 
    h1 = histogram(unitST, edges, 'Normalization','countdensity'); %Histogram of spike times binned in Xs bins, as specified by binw variable
    h1.FaceColor = [0.95 0.95 0.95]; 
end
%%
x = 13;
tmp = (mat(x).ccg_rest - mat(x).fr_rest)./mat(x).fr_rest;
figure;
for y = 1:size(tmp,2)
    subplot(5,3,y);
    plot(time, tmp(:,y));
    title(sprintf('%d %d',mat(x).m(y),mat(x).n(y)));
end