function ProjectData(pRFFolder, surfFolder, T1, Hemis)
% project volume data on original T1 surface using samsrf functions
%   Subject: subject ID
%   RootFolder: defaults to 'E:\Mojito\nifti\size\'
%   FSFolder: defaults to 'C:\Users\Ben\Desktop\FreeSurferSubjectFolder\'
% Rearranged for ERF study, 2019
% samsrf6.05, MATLAB R2016b, SPM12
% ----------------------------------------------------------------------- %
% prep
BatchNo=1;
% ----------------------------------------------------------------------- %
% let's do it
for iHemi = 1:length(Hemis)
    timeseries = dir([pRFFolder '*Avg_4D.nii']); % firstodd & theneven
    timeseries = {timeseries.name};
%     Structural=dir([mriFolder '*.nii']); %should be only one .nii in here
%     Structural=Structural.name;
    for i = 1:length(timeseries)
        timeserie = timeseries{i};
        cd(pRFFolder);
% Srf = samsrf_vol2srf(funimg,                strimg,  hemsurf,               ctxsteps, rule, nrmls)
        samsrf_vol2srf([pRFFolder timeserie], T1, [surfFolder Hemis{iHemi}], 0.5, '~', false, false); %here we don't normalise nor average because this was done when averaging niftis
    end%for iTS
end%for iHemi
end