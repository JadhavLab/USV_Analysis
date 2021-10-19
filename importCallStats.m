function callStats = importCallStats(workbookFile, sheetName, dataLines)
%IMPORTFILE Import data from a spreadsheet
%  C1P61BL20211007156PMSTATS = importCallStats(FILE) reads data from the
%  first worksheet in the Microsoft Excel spreadsheet file named FILE.
%  Returns the data as a table.
%
%  C1P61BL20211007156PMSTATS = importCallStats(FILE, SHEET) reads from the
%  specified worksheet.
%
%  C1P61BL20211007156PMSTATS = importCallStats(FILE, SHEET, DATALINES) reads
%  from the specified worksheet for the specified row interval(s).
%  Specify DATALINES as a positive scalar integer or a N-by-2 array of
%  positive scalar integers for dis-contiguous row intervals.
%
%  Example:
%  C1P61BL20211007156PMStats = importCallStats("G:\USV data\CallStats2021-10-15\C1_P6_1BL 2021-10-07  1_56 PM_Stats.xlsx", "Sheet1", [2, 485]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 15-Oct-2021 22:15:22

%% Input handling

% If no sheet is specified, read first sheet
if nargin == 1 || isempty(sheetName)
    sheetName = 1;
end

% If row start and end points are not specified, define defaults
if nargin <= 2
    dataLines = [2, 1000];
end

%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 17);

% Specify sheet and range
opts.Sheet = sheetName;
opts.DataRange = "A" + dataLines(1, 1) + ":Q" + dataLines(1, 2);

% Specify column names and types
opts.VariableNames = ["ID", "Label", "Accepted", "Score", "BeginTimes", "EndTimes", "CallLengths", "PrincipalFrequencykHz", "LowFreqkHz", "HighFreqkHz", "DeltaFreqkHz", "FrequencyStandardDeviationkHz", "SlopekHzs", "Sinuosity", "MeanPowerdBHz", "Tonality", "PeakFreqkHz"];
opts.VariableTypes = ["double", "string", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify variable properties
opts = setvaropts(opts, ["Label", "Accepted"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Label", "Accepted"], "EmptyFieldRule", "auto");

% Import the data
callStats = readtable(workbookFile, opts, "UseExcel", false);

for idx = 2:size(dataLines, 1)
    opts.DataRange = "A" + dataLines(idx, 1) + ":Q" + dataLines(idx, 2);
    tb = readtable(workbookFile, opts, "UseExcel", false);
    callStats = [callStats; tb]; %#ok<AGROW>
end

callStats=callStats(1:find(~isnan(callStats.ID),1,'last'),:); % remove trailing empty rows

end