%%
% N = X mice: DA to ACh pause
tmp = {}; for x = 1:length(beh); tmp{x} = strtok(beh(x).rec,'-'); end
uni = unique(tmp); nAn = length(uni);
da2peak_an = cell(nAn,3); ach2peak_an = cell(nAn,3);

r_amp = []; p_amp = []; r_lag = []; p_lag = [];
fit_amp = cell(1, 3); fit_lag = cell(1, 3);
slope_amp = []; slope_lag = [];

for y = 1:2
    for x = 1:nAn
        ii = find(strcmp(tmp,uni{x}));
        da2peak_an{x,y} = [da2achpeak{ii,y}]; % Concatenate all STAs
        ach2peak_an{x,y} = [ach2achpeak{ii,y}]; % Concatenate all STAs
        da2peak_an{x,y} = da2peak_an{x,y} - da2peak_an{x,y}(sta_time == -0.5,:);
        ach2peak_an{x,y} = ach2peak_an{x,y} - ach2peak_an{x,y}(sta_time == -0.5,:);
        
        max_ach = max(ach2peak_an{x,y});
        [mm_da1, ii] = max(da2peak_an{x,y}(find(sta_time == 0):find(sta_time == 0.18),:)); % find MAX within range of +0 to +0.3s after ACh peak maximal deflection
        max_time = sta_time(find(sta_time == 0) + ii - 1);
        mm_da2 = min(da2peak_an{x,y}(find(sta_time == -0.14):find(sta_time == 0),:)); % find MIN within range preceding ACh peak
        max_da = mm_da1 - mm_da2; % Difference in amplitude
        
        mdl_amp = fitlm(max_ach, max_da); % Linear regression: DA ~ ACh
        mdl_lag = fitlm(max_ach, max_time); 
        
        r_amp(x,y) = mdl_amp.Rsquared.Adjusted; r_lag(x,y) = mdl_lag.Rsquared.Adjusted;
        p_amp(x,y) = mdl_amp.Coefficients{2,4}; p_lag(x,y) = mdl_lag.Coefficients{2,4};
        fit_amp{y}(:,x) = mdl_amp.Coefficients{1,1} + mdl_amp.Coefficients{2,1}.*[0:30];
        fit_lag{y}(:,x) = mdl_lag.Coefficients{1,1} + mdl_lag.Coefficients{2,1}.*[0:30];
        slope_amp(x,y) = mdl_amp.Coefficients{2,1}; slope_lag(x,y) = mdl_lag.Coefficients{2,1};
    end
end

%% Single examples
fig = figure; fig.Position(3) = 1375;

subplot(1,3,1); hold on
% plot(max_ach, max_da, '.k');
% plot([0:30], mdl_amp.Coefficients{1,1} + mdl_amp.Coefficients{2,1}.*[0:30],'r','LineWidth',2)
plot([0:30], fit_amp{1}); 
xlabel('ACh peak amp (% dF/F)'); ylabel('DA peak amp (% dF/F)');
xlim([0 30]); ylim([-10 20]); axis('square');
title('DA peak amp ~ ACh peak amp (imm)');

subplot(1,3,2); hold on
% plot(max_ach, max_time, '.k');
% plot([0:30], mdl_lag.Coefficients{1,1} + mdl_lag.Coefficients{2,1}.*[0:30],'r','LineWidth',2)
plot([0:30], fit_lag{1}.*1000); 
xlabel('ACh peak amp (% dF/F)'); ylabel('DA peak latency (ms)'); ylim([0 200])
title('DA peak latency ~ ACh peak amp (imm)'); axis('square');

subplot(1,3,3);
violinplot([r_amp,r_lag]);
xticklabels({'amp imm','amp mov','lag imm','lag mov'});
ylabel('Rsquared'); ylim([-0.1 0.5]);
title('DA amp or lag ~ ACh peak amp'); axis('square');
movegui(gcf,'center');