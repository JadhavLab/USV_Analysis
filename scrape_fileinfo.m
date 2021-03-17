function [pd,animal_id] = scrape_fileinfo(filename)
    day = str2num(filename.name(1));
    pd = 4+2*day;
    
    %Ugly code to pull animal_id from file name
    %File name convention assumed to be e.g. 1_xf217-1BR Mar-16-2021  2_32 PM_Stats.xlsx
    first = strfind(filename.name,'-')+1;
    last = strfind(filename.name,' ')-1;
    fname = filename.name;
    animal_id = str2num(convertCharsToStrings(fname(first:last)));
end