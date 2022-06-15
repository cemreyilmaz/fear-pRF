function MergeData(SPMFolder, pRFFolder, NoOfRuns, FolderVector)
% merge volumes into 4D timeseries (one per run)
% Rearranged for ERF study, 2019
% samsrf6.05, MATLAB R2016b, SPM12
%   Subject: subject ID
%   NoOfRuns defaults to 6
%   RootFolder: defaults to 'E:\Mojito\nifti\'
%   FSFolder: defaults to 'C:\Users\Ben\Desktop\FreeSurferSubjectFolder\'
%   FolderVector: optional cellstr containing all subfolders except main exp runs
%   defaults to {'Hrf1', 'Scenes1', 'Scenes2', 'Rest1'}
% ----------------------------------------------------------------------- %
% prep
BatchNo=1;
for iRun=1:NoOfRuns
    FolderVector=[FolderVector, ['Ridge' num2str(iRun)]]; % append main exp runs
end
Subs=getsubfolders(SPMFolder); % should get ridge1, 2, 3 etc.
% ----------------------------------------------------------------------- %
% let's do it
for iFolder=1:length(FolderVector)
    Folder=FolderVector{iFolder};
    if  sum(strcmp(Subs, Folder))
        EPIs={};%initialise
        Folder=[SPMFolder Folder filesep];
        Files=dir([Folder '*.nii']);
        for iFile=1:length(Files)
            EPIs{iFile}=[Folder Files(iFile).name];
        end
        spm_file_merge(EPIs, [pRFFolder FolderVector{iFolder} '_4D']);
    end
end
end