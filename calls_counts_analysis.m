% info = inputdlg(["Start day","Interval","End day"], 'double');
% 
% start_day = str2num(cell2mat(info(1)));
% interval = str2num(cell2mat(info(2)));
% end_day = str2num(cell2mat(info(3)));

start_day = 6;  %First recordings on P6
interval = 2;   %Record every 2 days
end_day = 18;   %Last recordings on P18

%% initialize variables
m_duration = [];
m_delta_freq_kHz = [];
m_frequency_standard_deviation_kHz = [];
m_principal_frequency_kHz = [];
m_tonality = [];
m_slope_kHzs = [];
m_mean_power_dBHz = [];

label = [];

%% Get daily means for all variables

subj_means = [];

for k=start_day:interval:end_day
    dummy = (combined_usvs.Var1==int2str(k));
    subtable = combined_usvs(dummy,:);
    
    m_duration(end+1) = mean(subtable.CallLengths);
    
    
%     delta_freq_kHz(end+1) = mean(subtable.DeltaFreqkHz);
%     frequency_standard_deviation_kHz(end+1) = mean(subtable.FrequencyStandardDeviationkHz);
%     principal_frequency_kHz(end+1) = mean(subtable.PrincipalFrequencykHz);
%     tonality(end+1) = mean(subtable.Tonality);
%     slope_kHzs(end+1) = mean(subtable.SlopekHzs);
%     mean_power_dBHz(end+1) = mean(subtable.MeanPowerdBHz);
end

%boxplot(combined_usvs.CallLengths,combined_usvs.Var1)
%axis([0,8,0,.4]);

