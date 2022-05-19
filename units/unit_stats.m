fr_imm = nan(length(cinWT),1); 
cv_imm = nan(length(cinWT),1); 

for x = 1:length(cinWT)
    st = cinWT(x).st;
    ii = find(strcmp({behWT.rec},cinWT(x).rec));
    st = extractEventST(st, behWT(ii).onRest, behWT(ii).offRest, 0);
    isi = diff(st);
    fr_imm(x) = 1/mean(isi);
    cv_imm(x) = std(isi)/mean(isi);
end
cv_imm(114) = [];

%%
a = fr_imm;
fprintf('IMM firing rate: %1.2f +/- %1.2f Hz; range: %1.2f +/- %1.2f Hz; n = %d pCINs from %d recordings \n', nanmean(a), SEM(a,1), min(a), max(a), length(cinWT), length(unique({cinWT.rec})));

a = cv_imm;
fprintf('IMM coefficient of variation : %1.2f +/- %1.2f ; range: %1.2f +/- %1.2f \n', nanmean(a), SEM(a,1), min(a), max(a));

a = [cinWT.fr]';
fprintf('FULL firing rate: %1.2f +/- %1.2f Hz; range: %1.2f +/- %1.2f Hz\n', nanmean(a), SEM(a,1), min(a), max(a));

a = [cinWT.CV]'; a(114) = [];
fprintf('FULL coefficient of variation : %1.2f +/- %1.2f ; range: %1.2f +/- %1.2f\n', nanmean(a), SEM(a,1), min(a), max(a));