function [res] = readStationInfo(filename)
% READSTATIONINFO  Parses a RINEX station log file for receiver/antenna/monument change dates.
%
% INPUT:
%   filename : string, path to the site log file (e.g., 'TSKB00JPN_20221018.log')
%
% OUTPUT:
%   res      : struct with datetime arrays for fields:
%                  - rec : receiver change dates (datetime)
%                  - ant : antenna change dates (datetime)
%                  - mon : monument installation dates (datetime, if present)
%
% Each field contains a list of [yyyy mm dd] dates (converted to datetime).
%
% The function reads the file line-by-line, searching for lines with
% "Date Installed", and assigns them to the correct category according to the first character.
%
% Example:
%   info = readStationInfo('TSKB00JPN_20221018.log');
%   info.rec   % receiver change dates
%   info.ant   % antenna change dates

fid = fopen(filename);     % Open file for reading
counter = 1;              % Line counter
varAkt = '';              % Active variable (rec/ant/mon)
counters = struct();      % Counters for each type

while 1
    l = fgetl(fid);       % Read next line

    if l == -1
        % End of file reached
        warning(['file ', filename, ' ended after ', num2str(counter), ' lines!']);
        break
    end

    % Identify block type by first character in line
    if ~isempty(l)
        switch l(1)
            case '1'
                varAkt = 'mon';  % Monument
                if ~isfield(counters, varAkt)
                    counters.(varAkt) = 1;
                end
            case '3'
                varAkt = 'rec';  % Receiver
                if ~isfield(counters, varAkt)
                    counters.(varAkt) = 1;
                end
            case '4'
                varAkt = 'ant';  % Antenna
                if ~isfield(counters, varAkt)
                    counters.(varAkt) = 1;
                end
        end
    end

    % Check for "Date Installed" pattern in line
    if length(l) > 32
        switch l(1:31)
            case '     Date Installed           :'
                % Parse year, month, day from fixed columns
                yyyy = str2double(l(33:36));
                mm   = str2double(l(38:39));
                dd   = str2double(l(41:42));
                % Store date in appropriate field
                res.(varAkt)(counters.(varAkt),:) = [yyyy, mm, dd];
                counters.(varAkt) = counters.(varAkt) + 1;
        end
    end
    counter = counter + 1; % Increment line counter
end

fclose(fid); % Always close the file!

% Convert date arrays to datetime, removing rows with NaNs (if any)
if isfield(res, 'rec')
    res.rec = datetime(res.rec(~isnan(res.rec(:,1)),:));
end
if isfield(res, 'ant')
    res.ant = datetime(res.ant(~isnan(res.ant(:,1)),:));
end
if isfield(res, 'mon')
    res.mon = datetime(res.mon(~isnan(res.mon(:,1)),:));
end

end
