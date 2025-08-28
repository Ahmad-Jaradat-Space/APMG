function [configFileList] = getConfigsInFolder(varargin)
%[configFileList] = getConfigsInFolder(folderName) - Creates a cell with
%   the full paths to the configuration files. If no folder is given it is
%   assumed, that a subfolder named 'config' exists in the current
%   directory and that this is the target folder.
%
%   Input:      folderName...   Optional argument to read the contents of a
%                               specific folder
%
%   Output:     configFileList..Cell containing the full paths to the
%                               configuration files
%
%
%   @author: Matthias Gärtner
%   @date: 08.10.2019
%

narginchk(0,1)

if nargin<1
    folderName=[pwd, '/Config'];
else
    folderName = varargin{1};
end

dirContent = dir(fullfile(folderName, '*.txt'));

numberOfConfigs = length(dirContent);

configFileList=cell(numberOfConfigs,1);

for ii=1:numberOfConfigs
    configFileList{ii}=fullfile(folderName,dirContent(ii).name);    
end


end

