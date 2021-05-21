
%% Get basic info about testing days, number of subjects, etc.

info = inputdlg(["Start day","Interval","End day"], 'double');
 
start_day = str2num(cell2mat(info(1)));    %First recordings on P6
interval = str2num(cell2mat(info(2)));     %Record every 2 days
end_day = str2num(cell2mat(info(3)));      %Last recordings on P18
n_days = length(start_day:interval:end_day);

call_list={'count','CallLengths','PrincipalFrequencykHz',...
    'SlopekHzs','Tonality','DeltaFreqkHz'};
call_list=cohortfull.Properties.VariableNames;

[indx_params] = listdlg('ListString',call_list,...
    'PromptString',{'Select parameters to analyze'});

[indx_analysis] = listdlg('ListString',{'kruskal-wallis anova','mean ranks multiple comparisons test'...
    'plot subject means (error bars)','plot subject means (individuals)'},...
    'PromptString',{'Select analyses to run'});

[subj_names,a,b] = unique(cohortfull.ratNumber);
subj_geno=cohortfull.Genotype(a);
n_subjects = length(subj_names);
mycolormap = parula(length(subj_names));

% or genotype colormap
[a,b,allgeno]=unique(subj_geno);
mycolormap=lines(length(a));
mycolormap=mycolormap(allgeno,:);
% this goes over all

%% Counts analysis

for i=1:length(indx_params)

    subj_means=[];
    for k=1:n_days
        day = (start_day-interval)+(k*interval);
        indexes = cohortfull.Day==day;
        subT = cohortfull(indexes,:);               %Subtable containg only day k's data
        if indx_params(i) ==1
            
            for n=1:length(subj_names)
                indiv_indexes = subT.ratNumber==subj_names(n);
                subj_means(n,k) = sum(indiv_indexes);
                if sum(indiv_indexes)==0, subj_means(n,k)=nan; end
            end
           
        else
            for n=1:length(subj_names)
                indiv_indexes = subT.ratNumber==subj_names(n);
                subsubT = subT(indiv_indexes,:);
                subj_means(n,k) = mean(subsubT.(call_list{indx_params(i)}));
                subj_sds(n,k) = std(subsubT.(call_list{indx_params(i)}));
            end
        end
    end

    day_means=mean(subj_means);
    % Figures
    
    %K-W test
    if (any(indx_analysis==1))
        [p,tab,stats] = kruskalwallis(subj_means);
        hold on; box off;
        xticklabels([start_day:interval:end_day]);
        title('Number of calls');
        xlabel('Postnatal Day');
        ylabel('Number of USVs');
    end
    if (any(indx_analysis==2))
        multcompare(stats);
    end


    % Plot individual and overall means by day, with error bar
    if (any(indx_analysis==3))
        figure(1+i); hold on;
        errorbar(mean(subj_means,2,'omitnan'),SEM(subj_means,2),'k-','LineWidth',1)
            %plot(daily_subject_means,'k-','LineWidth',1)
        %plotSpread(day_counts);
        box off; zoom out;
        title('Number of calls');
        xlabel('Postnatal Day');
        ylabel('Number of USVs');
        xticklabels(start_day:interval:end_day);
    end

    % Plot individual means by day, tracking each subject
    if (any(indx_analysis==4))
        figure(30+i); hold on;
        day_counts_sorted = sortrows(subj_means,1,'ascend');

        for m=1:size(day_counts_sorted,1)
            pl(m)=plot(day_counts_sorted(m,:),'Color',mycolormap(m,:),'LineWidth',1.5);
        end
        box off; zoom out;
        title(call_list{indx_params(i)});
        xlabel('Postnatal Day');
        ylabel(call_list{indx_params(i)});
        xticklabels(start_day:interval:end_day);
        legend(pl(b),a);
    end
end

%% Duration analysis

if (any(indx_params==2))
    % initialize variables
    m_duration = [];
	se_duration = [];
    
    % Get daily means for all variables
    subj_means_duration = [];
    subj_sds = [];
    subj_names = unique(cohortfull.ratNumber);

    for k=1:n_days
        day = (start_day-interval)+(k*interval);
        indexes = cohortfull.day==int2str(day);
        subT = cohortfull(indexes,:);

        m_duration(end+1) = mean(subT.CallLengths);
        for n=1:length(subj_names)
            indiv_indexes = subT.ratNumber==subj_names(n);
            subsubT = subT(indiv_indexes,:);
            subj_means_duration(n,k) = mean(subsubT.CallLengths);
            subj_sds(n,k) = std(subsubT.CallLengths);
        end

        m_duration(end) = mean(subT.CallLengths); %Mean duration for all calls that day
    end

    se_duration = std(subj_sds, 'omitnan');
    daily_subject_means_duration = mean(subj_means_duration, 'omitnan');  %Mean across subjects per day
    subject_7day_duration = mean(subj_means_duration, 2, 'omitnan');      %Mean across days per subject
    
    % Figures
    
    %K-W test
    if (any(indx_analysis==1))
        [p,tab,stats] = kruskalwallis(subj_means_duration);
        hold on; box off;
        xticklabels([start_day:interval:end_day]);
        title('Call Duration');
        xlabel('Postnatal Day');
        ylabel('Mean subject call duration (sec)');
    end
    if (any(indx_analysis==2))
        multcompare(stats);
    end
        
    % Plot individual and overall means by day, with error bar
    if (any(indx_analysis==3))
        figure(20); hold on;
        errorbar(daily_subject_means_duration,se_duration,'k-','LineWidth',1)
        plotSpread(subj_means_duration)
        %plot(daily_subject_means,'k-','LineWidth',1)
        box off; zoom out;
        xticklabels([start_day:interval:end_day]);
        title('Call Duration');
        xlabel('Postnatal Day');
        ylabel('Mean subject call duration (sec)');
    end
        
    % Plot individual means by day, tracking each subject
    if (any(indx_analysis==4))
        figure(21); hold on;
        subj_means_duration_sorted = sortrows(subj_means_duration,1,'ascend');
        for m=1:length(subj_names)
            plot(subj_means_duration_sorted(m,:),'Color',mycolormap(m,:),'LineWidth',1.5);
        end
        title('Call Duration');
        xlabel('Postnatal Day');
        ylabel('Mean subject call duration (sec)');
        xticklabels([start_day:interval:end_day]);
    end
end

%% Frequency analysis

if (any(indx_params==3))
    % initialize variables
    m_principal = [];
    se_principal = [];
    m_delta = [];
    se_delta = [];
    
    % Get daily means for all variables
    subj_means_principal = [];
    subj_sds = [];
    subj_means_delta = [];
    subj_sds_delta = [];
    subj_names = unique(cohortfull.ratNumber);

    for k=1:n_days
        day = (start_day-interval)+(k*interval);
        indexes = cohortfull.day==int2str(day);
        subT = cohortfull(indexes,:);

        m_principal(end+1) = mean(subT.PrincipalFrequencykHz);
            m_delta(end+1) = mean(subT.DeltaFreqkHz);
        for n=1:length(subj_names)
            indiv_indexes = subT.ratNumber==subj_names(n);
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
    
    % Figures
    
    %K-W test
    if (any(indx_analysis==1))
        [p,tab,stats_pr] = kruskalwallis(subj_means_principal);
        hold on; box off;
        xticklabels([start_day:interval:end_day]);
        title('Principal USV frequency');
        xlabel('Postnatal Day');
        ylabel('Mean subject principal USV frequency (kHz)');
        
        [p,tab,stats_d] = kruskalwallis(subj_means_delta);
        hold on; box off;
        xticklabels([start_day:interval:end_day]);
        title('USV Delta frequency');
        xlabel('Postnatal Day');
        ylabel('Mean subject USV dFrequency (kHz)');
    end
    if (any(indx_analysis==2))
        multcompare(stats_pr);
        multcompare(stats_d);
    end
        
    % Plot individual and overall means by day, with error bar
    if (any(indx_analysis==3))
        figure(30); hold on;
        errorbar(daily_subject_means_principal,se_principal,'k-','LineWidth',1)
        plotSpread(subj_means_principal)
        %plot(daily_subject_means,'k-','LineWidth',1)
        box off; zoom out;
        xticklabels([start_day:interval:end_day]);
        title('Principal USV frequency');
        xlabel('Postnatal Day');
        ylabel('Mean subject principal USV frequency (kHz)');
        
        figure(31); hold on;
        errorbar(daily_subject_means_delta,se_delta,'k-','LineWidth',1)
        plotSpread(subj_means_delta)
        %plot(daily_subject_means,'k-','LineWidth',1)
        box off; zoom out;
        xticklabels([start_day:interval:end_day]);
        title('USV Delta frequency');
        xlabel('Postnatal Day');
        ylabel('Mean subject USV dFrequency (kHz)');
    end
        
    % Plot individual means by day, tracking each subject
    if (any(indx_analysis==4))
        figure(32); hold on;
        subj_means_principal_sorted = sortrows(subj_means_principal,1,'ascend');
        for m=1:12
            plot(subj_means_principal_sorted(m,:),'Color',colormap12(m,:),'LineWidth',1.5);
        end
        title('Principal USV frequency');
        xlabel('Postnatal Day');
        ylabel('Mean subject principal USV frequency (kHz)');
        xticklabels([start_day:interval:end_day]);

        figure(33); hold on;
        subj_means_delta_sorted = sortrows(subj_means_delta,1,'ascend');
        for m=1:length(subj_names)
            plot(subj_means_delta_sorted(m,:),'Color',mycolormap(m,:),'LineWidth',1.5);
        end
        title('USV Delta frequency');
        xlabel('Postnatal Day');
        ylabel('Mean subject USV dFrequency (kHz)');
        xticklabels([start_day:interval:end_day]);
    end
end

%% Slope analysis

if (any(indx_params==4))
    % initialize variables
    m_slope = [];
    se_slope = [];
    
    % Get daily means for all variables
    subj_means_slope = [];
    subj_sds = [];
    subj_names = unique(cohortfull.ratNumber);

    for k=1:n_days
        day = (start_day-interval)+(k*interval);
        indexes = cohortfull.day==int2str(day);
        subT = cohortfull(indexes,:);

        m_slope(end+1) = mean(subT.SlopekHzs);
        for n=1:length(subj_names)
            indiv_indexes = subT.ratNumber==subj_names(n);
            subsubT = subT(indiv_indexes,:);
            subj_means_slope(n,k) = mean(subsubT.SlopekHzs);
            subj_sds(n,k) = std(subsubT.SlopekHzs);
        end

        m_slope(end) = mean(subT.SlopekHzs); %Mean slope for all calls that day
    end

    se_slope = std(subj_sds, 'omitnan');
    daily_subject_means = mean(subj_means_slope, 'omitnan');  %Mean across subjects per day
    subject_7day = mean(subj_means_slope, 2, 'omitnan');      %Mean across days per subject
    
    % Figures
    
    %K-W test
    if (any(indx_analysis==1))
        [p,tab,stats] = kruskalwallis(subj_means_slope);
        hold on; box off;
        xticklabels([start_day:interval:end_day]);
        title('USV slope');
        xlabel('Postnatal Day');
        ylabel('Mean subject USV slope (kHz/s)');
    end
    if (any(indx_analysis==2))
        multcompare(stats);
    end    
        
    % Plot individual and overall means by day, with error bar
    if (any(indx_analysis==3))
        figure(40); hold on;
        errorbar(daily_subject_means,se_slope,'k-','LineWidth',1)
        plotSpread(subj_means_slope)
        %plot(daily_subject_means,'k-','LineWidth',1)
        box off; zoom out;
        title('USV slope');
        xlabel('Postnatal Day');
        ylabel('Mean subject USV slope (kHz/s)');
        xticklabels([start_day:interval:end_day]);
    end
        
    % Plot individual means by day, tracking each subject
    if (any(indx_analysis==4))
        figure(41); hold on;
        subj_means_slope_sorted = sortrows(subj_means_slope,1,'ascend');
        for m=1:length(subj_names)
            plot(subj_means_slope_sorted(m,:),'Color',mycolormap(m,:),'LineWidth',1.5);
        end
        title('USV slope');
        xlabel('Postnatal Day');
        ylabel('Mean subject USV slope (kHz/s)');
        xticklabels([start_day:interval:end_day]);
    end
end   
%% Tonality analysis

if (any(indx_params==5))
    % initialize variables
    m_tonality = [];
    se_tonality = [];

    % Get daily means for all variables
    subj_means_tonality = [];
    subj_sds = [];
    subj_names = unique(cohortfull.ratNumber);

    for k=1:n_days
        day = (start_day-interval)+(k*interval);
        indexes = cohortfull.day==int2str(day);
        subT = cohortfull(indexes,:);

        m_tonality(end+1) = mean(subT.Tonality);
        for n=1:length(subj_names)
            indiv_indexes = subT.ratNumber==subj_names(n);
            subsubT = subT(indiv_indexes,:);
            subj_means_tonality(n,k) = mean(subsubT.Tonality);
            subj_sds(n,k) = std(subsubT.Tonality);
        end

        m_tonality(end) = mean(subT.Tonality); %Mean tonality for all calls that day
    end

    se_tonality = std(subj_sds, 'omitnan');
    daily_subject_means = mean(subj_means_tonality, 'omitnan');  %Mean across subjects per day
    subject_7day = mean(subj_means_tonality, 2, 'omitnan');      %Mean across days per subject
    
    % Figures
    
    %K-W test
    if (any(indx_analysis==1))
        [p,tab,stats] = kruskalwallis(subj_means_tonality);
    hold on; box off;
    xticklabels([start_day:interval:end_day]);
    title('Tonality');
    xlabel('Postnatal Day');
    ylabel('Mean subject USV tonality (dB/kHz)');
    end
    if (any(indx_analysis==2))
        multcompare(stats);
    end
        
    % Plot individual and overall means by day, with error bar
    if (any(indx_analysis==3))
        figure(50); hold on;
        errorbar(daily_subject_means,se_tonality,'k-','LineWidth',1)
        plotSpread(subj_means_tonality)
        %plot(daily_subject_means,'k-','LineWidth',1)
        box off; zoom out;
        title('Tonality');
        xlabel('Postnatal Day');
        ylabel('Mean subject USV tonality (dB/kHz)');
        xticklabels([start_day:interval:end_day]);
    end
        
    % Plot individual means by day, tracking each subject
    if (any(indx_analysis==4))
        figure(51); hold on;
        subj_means_tonality_sorted = sortrows(subj_means_tonality,1,'ascend');
        for m=1:length(subj_names)
            plot(subj_means_tonality_sorted(m,:),'Color',mycolormap(m,:),'LineWidth',1.5);
        end
        title('Tonality');
        xlabel('Postnatal Day');
        ylabel('Mean subject USV tonality (dB/kHz)');
        xticklabels([start_day:interval:end_day]);
    end
end