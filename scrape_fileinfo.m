%% Scrape file info

% Pulls animal ID and day tested from file names
% File name convention assumed to be e.g. 1_xf217-1BR Mar-16-2021  2_32 PM_Stats.xlsx
% Modify as needed based on your file name template
%P8_2RR May-19-2021  5_37 PM_Stats.xlsx


function [pd,animal_id,fname] = scrape_fileinfo(filename)
    day = str2num(filename.name(1));
    pd = 4+2*day;
   
    %Ugly code to pull animal_id from file name
    first = strfind(filename.name,'-')+1;
    last = strfind(filename.name,' ')-1;
    fname = filename.name;
    animal_id = str2num(convertCharsToStrings(fname(first(1):last(1))));
    

end