function [dahiti_combined, theia_combined] = read_and_combine_time_series(base_directory)
    % Define paths for Dahiti and Theia data
    dahiti_path = fullfile(base_directory, 'dahiti', 'Rhine_River');
    theia_path = fullfile(base_directory, 'thalia', 'Rhine_River');
    
    % Initialize empty tables to store combined time series
    dahiti_combined = table();
    theia_combined = table();
    
    % Read Dahiti data (txt files in subfolders)
    dahiti_folders = dir(fullfile(dahiti_path, '*'));
    dahiti_folders = dahiti_folders([dahiti_folders.isdir]); % Only keep folders
    dahiti_folders = dahiti_folders(~ismember({dahiti_folders.name}, {'.', '..'})); % Exclude . and ..
    
    for i = 1:length(dahiti_folders)
        folder = fullfile(dahiti_path, dahiti_folders(i).name);
        txt_files = dir(fullfile(folder, '*.txt'));
        
        for j = 1:length(txt_files)
            file_path = fullfile(folder, txt_files(j).name);
            
            % Read Dahiti data
            fid = fopen(file_path, 'r');
            data = textscan(fid, '%s %f %f', 'CommentStyle', '#', 'HeaderLines', 10); % Skip rows starting with #
            fclose(fid);
            
            % Convert to table
            data_table = table(data{1}, data{2}, data{3}, 'VariableNames', {'Date', 'WaterLevel', 'Error'});
            data_table.Date = datetime(data_table.Date, 'InputFormat', 'yyyy-MM-dd'); % Convert to datetime
            dahiti_combined = [dahiti_combined; data_table];
        end
    end
    
    % Read Theia data (csv files directly in the folder)
    theia_files = dir(fullfile(theia_path, '*.csv'));
    
    for i = 1:length(theia_files)
        file_path = fullfile(theia_path, theia_files(i).name);
        
        % Read Theia data
        fid = fopen(file_path, 'r');
        data = textscan(fid, '%s %*s %f %f %*[^\n]', 'CommentStyle', '#', 'HeaderLines', 31); % Skip rows starting with #
        fclose(fid);
        
        % Convert to table
        data_table = table(data{1}, data{2}, data{3}, 'VariableNames', {'Date', 'WaterLevel', 'Error'});
        data_table.Date = datetime(data_table.Date, 'InputFormat', 'yyyy-MM-dd'); % Convert to datetime
        theia_combined = [theia_combined; data_table];
    end
    
    % Sort the combined data by date
    dahiti_combined = sortrows(dahiti_combined, 'Date');
    theia_combined = sortrows(theia_combined, 'Date');
    
    % Display the first few rows of the combined data
    disp('Dahiti Combined Data:');
    disp(head(dahiti_combined));
    disp('Theia Combined Data:');
    disp(head(theia_combined));
end