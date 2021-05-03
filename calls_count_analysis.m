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
m_count = [];
se_count = [];
daily_means = [];
day_counts = [];

%% Get daily means for all variables

subj_names = unique(T.Var2);

for k=1:n_days
    day = (start_day-interval)+(k*interval);
    indexes = T.Var1==int2str(day);
    subT = T(indexes,:);               %Subtable containg only day k's data
    
    for n=1:length(subj_names)
        indiv_indexes = subT.Var2==subj_names(n);
        day_counts(n,k) = sum(indiv_indexes==1);
    end
end

daily_means = mean(day_counts, 'omitnan');      %Mean across subjects per day
subject_7day = mean(day_counts, 2, 'omitnan');  %Mean across days per subject
se_count = std(day_counts,1,'omitnan');

%% Figures


[p,tab,stats] = kruskalwallis(day_counts);
hold on;
box off;
xticklabels([start_day:interval:end_day]);
title('Number of calls');
xlabel('Postnatal Day');
ylabel('Number of USVs');

%% Plots

% Plot individual and overall means by day, with error bar
figure(5);
hold on;
errorbar(daily_means,se_count,'k-','LineWidth',1)
    %plot(daily_subject_means,'k-','LineWidth',1)
box off;
zoom out;
xticklabels([start_day:interval:end_day]);
title('Number of calls');
xlabel('Postnatal Day');
ylabel('Number of USVs');

% Plot individual means by day, tracking each subject
figure(8)
hold on;
colormap = parula;
colormap12 = colormap(1:256/12:256,:)
day_counts_sorted = sortrows(day_counts,1,'ascend')

for m=1:12
    plot(day_counts_sorted(m,:),'Color',colormap12(m,:),'LineWidth',1.5);
end
title('Number of calls');
xlabel('Postnatal Day');
ylabel('Number of USVs');
xticklabels([start_day:interval:end_day]);
