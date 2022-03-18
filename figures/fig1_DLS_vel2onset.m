a = [1:length(modACh)]; rmv = [3 17 20 22 23 30 36:39 41 48 50]; a(rmv) = [];
beh = modACh(a); % Extract recordings with reward

%%
align_vel = cell(length(beh),1);
Fs = beh(1).Fs; 
for x = 1:length(beh)
    ev = nan(length(beh(x).on),1);
    if all(logical(~rem(beh(x).on,1))); diffFs = 1; else; diffFs = 50; end
%     for ii = 1:length(beh(x).on) % iterate over movement bouts
%         acc = getAcc(beh(x).vel(beh(x).on(ii)*diffFs:beh(x).off(ii)*diffFs)); % extract acceleration during specified movement bout
%         [~,locs] = findpeaks(acc,'MinPeakProminence',0.5,'MinPeakDistance',0.5); % location of peaks, using findpeaks function
%         if ~isempty(locs)
%             locs_sec = beh(x).time(locs + beh(x).on(ii)*diffFs - 1); % acceleration peaks during movement bout, in seconds
%             ev(ii) = locs_sec(1); % first acceleration in movement bout, in seconds
%         end
%     end
    ev = beh(x).on./(Fs/diffFs);
    sig = beh(x).vel; % Signal that will be aligned to event times
    [mat,vel_time] = getSTA(sig, ev, Fs, [-2 1]);
    align_vel{x} = mat;
end

vel_adj = align_vel;
for x = 1:length(beh); y = 1;
    vel_adj{x,y} = vel_adj{x,y} - nanmean(vel_adj{x,y}(find(vel_time == -2):find(vel_time == -1),:));
end

% n = X mice
tmp = {}; for x = 1:length(beh); tmp{x} = strtok(beh(x).rec,'-'); end
uni = unique(tmp); nAn = length(uni);
vel_an = cell(length(uni),1); vel_an_avg = [];
for x = 1:nAn
    ii = find(strcmp(tmp,uni{x}));
    vel_an_avg(:,x) = nanmean([vel_adj{ii}],2);
    vel_an{x} = [vel_adj{ii}];
end

%% Population AVERAGE
figure;
shadederrbar(vel_time, nanmean(vel_an_avg,2), SEM(vel_an_avg,2), 'k');
xlabel('Latency to 1st Acceleration (s)'); xlim([-1 1]); xticks([-1:0.5:1]);
ylabel('Velocity (cm/s)'); % ylim([-3 7]); yticks([-2:2:6]);
title(sprintf('ACh3.0 DLS Population Average (n = %d mice)',size(align_an_avg,2)));