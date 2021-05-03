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
m_slope = [];
se_slope = [];

%% Get daily means for all variables

subj_means_slope = [];
subj_sds = [];
subj_names = unique(T.Var2);

for k=1:n_days
    day = (start_day-interval)+(k*interval);
    indexes = T.Var1==int2str(day);
    subT = T(indexes,:);
    
    m_slope(end+1) = mean(subT.SlopekHzs);
    for n=1:length(subj_names)
        indiv_indexes = subT.Var2==subj_names(n);
        subsubT = subT(indiv_indexes,:);
        subj_means_slope(n,k) = mean(subsubT.SlopekHzs);
        subj_sds(n,k) = std(subsubT.SlopekHzs);
    end
    
    m_slope(end) = mean(subT.SlopekHzs); %Mean slope for all calls that day
end

se_slope = std(subj_sds, 'omitnan');
daily_subject_means = mean(subj_means_slope, 'omitnan');  %Mean across subjects per day
subject_7day = mean(subj_means_slope, 2, 'omitnan');      %Mean across days per subject
%% Figures

[p,tab,stats] = kruskalwallis(subj_means_slope);
hold on; box off;
xticklabels([start_day:interval:end_day]);
title('USV slope');
xlabel('Postnatal Day');
ylabel('Mean subject USV slope (kHz/s)');

figure(3); hold on;
errorbar(daily_subject_means,se_slope,'k-','LineWidth',1)
plotSpread(subj_means_slope)
%plot(daily_subject_means,'k-','LineWidth',1)
box off; zoom out;
title('USV slope');
xlabel('Postnatal Day');
ylabel('Mean subject USV slope (kHz/s)');
xticklabels([start_day:interval:end_day]);


figure(7)
hold on;
colormap = parula;
colormap12 = colormap(1:256/12:256,:)

