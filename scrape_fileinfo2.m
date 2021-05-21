%% Scrape file info

% Pulls animal ID and day tested from file names
% File name convention assumed to be e.g. 1_xf217-1BR Mar-16-2021  2_32 PM_Stats.xlsx
% Modify as needed based on your file name template
%P8_2RR May-19-2021  5_37 PM_Stats.xlsx


function [pd,animal_id,filename] = scrape_fileinfo2(filestruct)
    filename=filestruct.name;
    split=find(filename=='_',1,'first');
    pd=filename(2:split-1);
    endcond=find(filename==' ',1,'first');
    animal_id = filename(split+1:endcond-1);
end