function [dirList] = getsubfolders(dir_name)

    %get_directory_names; this function outputs a cell with directory names (as
    %strings), given a certain dir name (string)
    %from: http://stackoverflow.com/questions/8748976/list-the-subfolders-
    %in-a-folder-matlab-only-subfolders-not-files
    % samsrf6.05, MATLAB R2016b, SPM12

    dd = dir(dir_name);
    isub = [dd(:).isdir]; % returns logical vector
    dirList = {dd(isub).name}';
    dirList(ismember(dirList,{'.','..'})) = [];

end