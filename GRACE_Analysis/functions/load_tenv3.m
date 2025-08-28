function [out] = load_tenv3(filename)

% load tenv3 files and extract data

%% Initialize variables.
startRow = 2;
endRow = inf;

% For more information, see the TEXTSCAN documentation.
formatSpec = '%4C%8s%10f%6f%5f%2f%7f%7f%10f%10f%10f%6f%10f%8f%9f%9f%9f%10f%10f%10f%15f%16f%12f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% Here, iterations over all lines in the given file are performed, within
% each of those all columns are extracted in 'dataArrayBlock'.
% Then, the second 'for' loop fills the lines of the cell array
% 'dataArray', from which we can extract any columns in the next part.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Remove white space around all cell columns.
dataArray{2} = strtrim(dataArray{2});

%% Close the text file.
fclose(fileID);

%% Allocate imported array to column variable names
% Extract the data from relevant colums. Refer to
% http://geodesy.unr.edu/gps_timeseries/README_tenv3load.txt and
% http://geodesy.unr.edu/gps_timeseries/README_tenv3.txt for the column
% content.

% tp = time vector in decimal year format
out.tp = dataArray{:,3};
% t = MJD
out.t = dataArray{:,4};
% up, east, north = coordinate time series
out.up = dataArray{:, 12} + dataArray{:, 13};
out.east = dataArray{:, 8} + dataArray{:, 9};
out.north = dataArray{:, 10} + dataArray{:, 11};
% formal error
out.sig_em = dataArray{:, 15};
out.sig_nm = dataArray{:, 16};
out.sig_um = dataArray{:, 17};

end