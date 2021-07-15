%% Scrape file info

% Pulls animal ID and day tested from file names
% File name convention assumed to be e.g. 1_xf217-1BR Mar-16-2021  2_32 PM_Stats.xlsx
% Modify as needed based on your file name template
%P8_2RR May-19-2021  5_37 PM_Stats.xlsx


function [cohort,animal_id,age] = scrape_fileinfo2(filename)
    % first trim to first space
    filename=filename(1:find(filename==' ',1,'first'));
    split=find(filename=='_');
    if length(split)~=2
        split=find(filename=='-');
    end
    cohort= filename(1:split(1)-1);
    age=filename(split(1)+1:split(2)-1);
    endcond=find(filename==' ',1,'first');
    animal_id = filename(split(2)+1:endcond-1);
end