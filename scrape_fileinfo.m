function [pd,animal_id] = scrape_fileinfo(filename)
    day = str2num(filename.name(1));
    pd = 4+2*day;
    
    %Ugly code to pull animal_id from file name
    %File name convention assumed to be e.g. 1_xf217-1BR Mar-16-2021  2_32 PM_Stats.xlsx
    dashdelim = strfind(filename.name,'-');
    spacedelim = strfind(filename.name,' ');
    fname = filename.name;
    animal_id = fname(dashdelim(1)+1:spacedelim(1)-1);
    animal_id = convertCharsToStrings(animal_id);
end