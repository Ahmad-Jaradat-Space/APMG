function config=readConfig(filename)
%   last change: 10.Apr 2018

fid=fopen(filename);

line=fgetl(fid);

% valid chars
va_char='[a-zA-Z0-9_.+-\/]';
config = struct;

%Schleife über alle Zeilen
while ischar(line)
    %ignore empty lines
    if ~isempty(line)
        % ignore comments %
        if ~strcmp(line(1),'%')
            tokens = regexp( line,['^\s*([a-zA-z0-9_]+)\s*=\s*(',va_char,'+(?:\s*,\s*',va_char,'+)*)\s*(?:%.*)?$'],'tokens');
            if isempty(tokens)
                warning('The line "%s" does not  match',line)
            else
                values = regexp(tokens{1}{2},'\s*\s','split');
                
                if size(values,2)>1
                    config.(tokens{1}{1})=values;
                else
                    config.(tokens{1}{1})=values{1};
                end
            end
            
            
        end
    end
    line=fgetl(fid);
end
fclose(fid);
end
