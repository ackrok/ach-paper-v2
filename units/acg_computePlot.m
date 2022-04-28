%% Compute ACG
% st = {data.clusters(data.lbl.cin).spikeTimes}; %extract spike times from data structure
st = {sub.st};

acgLin = []; %initialize matrixes
for x = 1:length(st)
    [xLin, nLin, ~, ~]= myACG(st{x},[],[]); %compute autocorrelogram with linear bins
    nLin = [flip(nLin);nLin]; %symmetrical ACG
    nLin = nLin./length(st{x});
    xLin = [-1.*flip(xLin);xLin]; xLin = xLin.*1000; %x-values for ACG in milliseconds
    acgLin = [acgLin, nLin]; %add to matrix
end

%% PLOT ACGs in subplots
figure;
plM = floor(sqrt(length(st))); plN = ceil(length(st)/plM); c = hsv(length(st)); 
for x = 1:length(st)
    subplot(plM,plN,x); hold on
    bar(xLin, acgLin(:,x),'FaceColor',c(x,:),'FaceAlpha',0.5);
    xlim([-250 250]); xlabel('Lag (ms)');
    title(sprintf('%d',x)); 
    % title(sprintf('ACG unit#%d',data.lbl.cin(x))); 
end

%% Plot overlapping ACGs
figure; hold on
c = hsv(length(st)); 
for x = 1:length(st)
    plot(xLin, acgLin(:,x),'Color',c(x,:),'DisplayName',sprintf('unit#%d',data.lbl.cin(x)));
end
xlim([-250 250]); xlabel('Lag (ms)');
title(sprintf('%s-%s: ACG',data.mouse,data.date)); 

%% Plot all CCGs from mat structuer
mat = struct; %NAME BASED ON COMB MAT STRUCTURE

figure; 
plM = floor(sqrt(length(mat))); plN = ceil(length(mat)/plM); c = hsv(length(mat)); 
for x = 1:length(mat)
    [xLin, nLin, ~, ~]= myACG(mat(x).st,[],[]); %compute autocorrelogram
    nLin = [flip(nLin);nLin]; %symmetrical ACG
    xLin = [-1.*flip(xLin);xLin]; xLin = xLin.*1000; %x-values for ACG in milliseconds
    
    subplot(plM,plN,x); 
    bar(xLin, nLin, 'FaceColor',c(x,:), 'FaceAlpha',0.5);
    xlim([-250 250]); xlabel('Lag (ms)');
end

%% SHADEDERRBAR
figure;
sm = 10; ds = 1;
shadederrbar(xLin(1:ds:end), movmean(nanmean(acgLin(1:ds:end,:),2),sm), movmean(SEM(acgLin(1:ds:end,:),2),sm), 'b');
ylabel('Firing Rate (Hz)');
xlabel('Lag (ms)'); xlim([-250 250]); 

%% complex shadederrbar
hold on;
sm = 10; ds = 1; color = 'k';

x = xLin(1:ds:end);
y = movmean(nanmean(acgLin(1:ds:end,:),2),sm); y = y - min(y);
z = movmean(SEM(acgLin(1:ds:end,:),2),sm);
if(size(x,1)~=1); x = x'; end
if(size(y,1)~=1); y = y'; end
if(size(z,1)~=1); z = z'; end
color=char2rgb(color);

patchcolor=color+(1-color)*.3; %Creates the patch color
yerru=y+z;
yerrl=y-z;
xpatch=[x,fliplr(x)]; %Creates x axis for the path
ypatch=[yerru,fliplr(yerrl)]; %Creates y axis for patch
hold on
fill(xpatch,ypatch,patchcolor,'FaceAlpha',0.5,'EdgeAlpha',0,'EdgeColor',patchcolor,'LineStyle','none'); %Creates the patch

patchcolor=color+(1-color)*.8; %Creates the patch color
yerru=y-z;
yerrl=zeros(1,length(x));
xpatch=[x,fliplr(x)]; %Creates x axis for the path
ypatch=[yerru,fliplr(yerrl)]; %Creates y axis for patch
hold on
fill(xpatch,ypatch,patchcolor,'FaceAlpha',0.5,'EdgeAlpha',0,'EdgeColor',patchcolor,'LineStyle','none'); %Creates the patch
main = plot(x,y,'-','Color',color); %Plots main data   