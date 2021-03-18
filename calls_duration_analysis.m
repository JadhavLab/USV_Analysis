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
m_duration = [];
se_duration = [];

%% Get daily means for all variables

subj_means = [];
subj_sds = [];
subj_names = unique(T.Var2);

for k=1:n_days
    day = (start_day-interval)+(k*interval);
    indexes = T.Var1==int2str(day);
    subT = T(indexes,:);
    
    m_duration(end+1) = mean(subT.CallLengths);
    for n=1:length(subj_names)
        indiv_indexes = subT.Var2==subj_names(n);
        subsubT = subT(indiv_indexes,:);
        subj_means(n,k) = mean(subsubT.CallLengths);
        subj_sds(n,k) = std(subsubT.CallLengths);
    end
    
    m_duration(end) = mean(subT.CallLengths); %Mean duration for all calls that day
end

se_duration = std(subj_sds, 'omitnan');
daily_subject_means = mean(subj_means, 'omitnan');

%% Figures

errorbar(daily_subject_means,se_duration,'k-','LineWidth',1)
hold on;
box off;
zoom out;
xticklabels([start_day:interval:end_day]);
title('Call Duration');
xlabel('Postnatal Day');
ylabel('Mean subject call duration (sec)');

p = anova1(subj_means);
hold on;
box off;
xticklabels([start_day:interval:end_day]);
title('Call Duration');
xlabel('Postnatal Day');
ylabel('Mean subject call duration (sec)');


%% daily duration means by animal
dayratduration=[];
for k=1:n_days
    day = (start_day-interval)+(k*interval);
    dayidx=T.Var1==int2str(day);
    [~,~,ratIDs]=unique(T.Var2);
    dayduration(:,k)=accumarray(ratIDs(dayidx), T.CallLengths(dayidx),[], @mean);
end

figure; plot(start_day:interval:end_day,dayduration','Color',[.7 .7 .7]);
hold on;
errorbar(start_day:interval:end_day,nanmean(dayduration),nanstd(dayduration),...
    'CapSize',0,'Color','k','LineWidth',4);
xlabel('Postnatal day'); ylabel('Mean Call Duration');

%% Howabout DeltaFreq

daydeltafreq=[];
for k=1:n_days
    day = (start_day-interval)+(k*interval);
    dayidx=T.Var1==int2str(day);
    [~,~,ratIDs]=unique(T.Var2);
    daydeltafreq(:,k)=accumarray(ratIDs(dayidx), T.DeltaFreqkHz(dayidx),[], @mean);
end

figure; plot(start_day:interval:end_day,daydeltafreq','Color',[.7 .7 .7]);
hold on;
errorbar(start_day:interval:end_day,nanmean(daydeltafreq),nanstd(daydeltafreq),...
    'CapSize',0,'Color','k','LineWidth',4);
xlabel('Postnatal day'); ylabel('Call Frequency Span');

%% and call count

daydeltafreq=[];
for k=1:n_days
    day = (start_day-interval)+(k*interval);
    dayidx=T.Var1==int2str(day);
    [~,~,ratIDs]=unique(T.Var2);
    daydeltafreq(:,k)=accumarray(ratIDs(dayidx), T.DeltaFreqkHz(dayidx),[], @length);
end

figure; plot(start_day:interval:end_day,daydeltafreq','Color',[.7 .7 .7]);
hold on;
errorbar(start_day:interval:end_day,nanmean(daydeltafreq),nanstd(daydeltafreq),...
    'CapSize',0,'Color','k','LineWidth',4);
xlabel('Postnatal day'); ylabel('Call count');

