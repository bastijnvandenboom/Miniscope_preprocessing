% normcorre_choose_data
% select file for NormCorre - non rigid movement correction
% 
% Bastijn van den Boom

%% select file 
if ~exist('name', 'var') | isempty(name) | name == 0
    % choose files manually 
    
    dir_nm = [cd(), filesep]; %use the current path
    [file_nm, dir_nm] = uigetfile(fullfile(dir_nm, '*.tif;*.mat;*.h5'));
    name = [dir_nm, file_nm];  % full name of the data file

    if isempty(name) | name == 0
        fprintf('No file was selected. Try again!\n');
    else
        fprintf('Successfully selected file %s\nDirectory: %s\n', file_nm, name);
    end
    
else
    % use pre-specified file 
    if exist('name', 'file')
        [dir_nm, file_nm, file_type] = fileparts(name);
    else
        dir_nm = 0; 
    end
end

