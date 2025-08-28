function dataStruct = read_all_time_series(baseDir, datasetType)
    % READ_ALL_TIME_SERIES Reads all Theia (CSV) or Dahiti (TXT) time-series in a directory
    % INPUT:
    %    baseDir - Base directory (e.g., 'Rhine_river/')
    %    datasetType - 'Theia' or 'Dahiti'
    % OUTPUT:
    %    dataStruct - Structure with fields for each station, containing time-series data

    % Initialize structure to store all time-series data
    dataStruct = struct();

    % Set file extension based on dataset type
    if strcmp(datasetType, 'Theia')
        fileExt = '*.csv';
        fileList = dir(fullfile(baseDir, fileExt));  % List all CSV files in base directory
    elseif strcmp(datasetType, 'Dahiti')
        fileExt = '*.txt';
        fileList = dir(fullfile(baseDir, '**', fileExt)); % List all TXT files in subdirectories
    else
        error('Unsupported dataset type. Choose ''Theia'' or ''Dahiti''.');
    end

    % Loop through each file
    for i = 1:length(fileList)
        filePath = fullfile(fileList(i).folder, fileList(i).name);
        stationName = erase(fileList(i).name, fileExt(2:end));  % Remove extension

        % Read time-series data
        data = read_time_series(filePath);

        % Store in structure with station name as field
        dataStruct.(stationName) = data;
    end

    disp(['Loaded ', num2str(length(fileList)), ' ', datasetType, ' time-series files.']);
end
