% info = inputdlg(["Start day","Interval","End day"], 'double');
% 
% start_day = str2num(cell2mat(info(1)));
% interval = str2num(cell2mat(info(2)));
% end_day = str2num(cell2mat(info(3)));

start_day = 6;  %First recordings on P6
interval = 2;   %Record every 2 days
end_day = 18;   %Last recordings on P18
n_days = length(start_day:interval:end_day);

%% initialize variables
m_principal = [];
se_principal = [];
m_delta = [];
se_delta = [];

%% Get daily means for all variables

subj_means_principal = [];
subj_sds = [];
subj_means_delta = [];
subj_sds_delta = [];
subj_names = unique(T.Var2);

for k=1:n_days
    day = (start_day-interval)+(k*interval);
    indexes = T.Var1==int2str(day);
    subT = T(indexes,:);
    
    m_principal(end+1) = mean(subT.PrincipalFrequencykHz);
        m_delta(end+1) = mean(subT.DeltaFreqkHz);
    for n=1:length(subj_names)
        indiv_indexes = subT.Var2==subj_names(n);
        subsubT = subT(indiv_indexes,:);
        subj_means_principal(n,k) = mean(subsubT.PrincipalFrequencykHz);
        subj_sds(n,k) = std(subsubT.PrincipalFrequencykHz);
        subj_means_delta(n,k) = mean(subsubT.DeltaFreqkHz);
        subj_sds_delta(n,k) = std(subsubT.DeltaFreqkHz);
        
    end  
    m_principal(end) = mean(subT.PrincipalFrequencykHz); %Mean principal freq for all calls that day
    m_delta(end) = mean(subT.DeltaFreqkHz); %Mean dFreq for all calls that day

end

se_principal = std(subj_sds, 'omitnan');
se_delta = std(subj_sds, 'omitnan');
daily_subject_means_principal = mean(subj_means_principal, 'omitnan');  %Mean across subjects per day
subject_7day_principal = mean(subj_means_principal, 2, 'omitnan');      %Mean across days per subject

daily_subject_means_delta = mean(subj_means_delta, 'omitnan');  %Mean across subjects per day
subject_7day_delta = mean(subj_means_delta, 2, 'omitnan');      %Mean across days per subject

%% Figures for principal frequency

[p,tab,stats] = kruskalwallis(subj_means_principal);
hold on; box off;
xticklabels([start_day:interval:end_day]);
title('Principal USV frequency');
xlabel('Postnatal Day');
ylabel('Mean subject principal USV frequency (kHz)');

figure(3); hold on;
errorbar(daily_subject_means_principal,se_principal,'k-','LineWidth',1)
plotSpread(subj_means_principal)
%plot(daily_subject_means,'k-','LineWidth',1)
box off; zoom out;
xticklabels([start_day:interval:end_day]);
title('Principal USV frequency');
xlabel('Postnatal Day');
ylabel('Mean subject principal USV frequency (kHz)');

%% Figures for delta frequency

[p,tab,stats] = kruskalwallis(subj_means_delta);
hold on; box off;
xticklabels([start_day:interval:end_day]);
title('USV Delta frequency');
xlabel('Postnatal Day');
ylabel('Mean subject USV dFrequency (kHz)');

figure(6);hold on;
errorbar(daily_subject_means_delta,se_delta,'k-','LineWidth',1)
plotSpread(subj_means_delta)
%plot(daily_subject_means,'k-','LineWidth',1)
box off; zoom out;
xticklabels([start_day:interval:end_day]);
title('USV Delta frequency');
xlabel('Postnatal Day');
ylabel('Mean subject USV dFrequency (kHz)');

figure(6)
hold on;
ylabel('Mean subject USV dFrequency (kHz)');

figure(9)
hold on;
colormap = parula;
colormap12 = colormap(1:256/12:256,:)
subj_means_principal_sorted = sortrows(subj_means_principal,1,'ascend')

for m=1:12
    plot(subj_means_principal_sorted(m,:),'Color',colormap12(m,:),'LineWidth',1.5);
end
title('Principal USV frequency');
xlabel('Postnatal Day');
ylabel('Mean subject principal USV frequency (kHz)');
xticklabels([start_day:interval:end_day]);

figure(10)
hold on;
colormap = parula;
colormap12 = colormap(1:256/12:256,:)
subj_means_delta_sorted = sortrows(subj_means_delta,1,'ascend')

for m=1:12
    plot(subj_means_delta_sorted(m,:),'Color',colormap12(m,:),'LineWidth',1.5);
end
title('USV Delta frequency');
xlabel('Postnatal Day');
ylabel('Mean subject USV dFrequency (kHz)');
xticklabels([start_day:interval:end_day]);