function pRF_DoSubject_Ridge(pRFFolder, FSFolder, surfFolder, SrfFiles, Rois)
% Rearranged for ERF study, 2019
% samsrf6.05, MATLAB R2016b, SPM12
cd (pRFFolder);
MakeOccRoi(surfFolder);
% ----------------------------------------------------------------------- %
% Check that everything is there
if ~exist(fullfile(cd, 'src_pRF.mat'), 'file')
    disp('Search space doesn''t exist.');
end
if ~exist(fullfile(cd, 'aps_Ridge.mat'), 'file')
    disp('Aps file doesn''t exist.');
    copyfile([FSFolder filesep 'aps_Ridge.mat'], ...
        [pRFFolder 'aps_Ridge.mat']);
    disp('Copied aps file from freesurfer directory');
end
%if ~exist(fullfile(cd, 'lh_occ.label'), 'file')
%    disp('Left hemisphere ROI label doesn''t exist!');
%end
%if ~exist(fullfile(cd, 'rh_occ.label'), 'file')
%    disp('Right hemisphere ROI label doesn''t exist!');
%end
% ----------------------------------------------------------------------- %
% Fit model for each hemisphere
for iHemi = 1:2
    Ridge_2D_Gaussian_Prf(SrfFiles(:,iHemi), pRFFolder, Rois{iHemi})
    end
end
