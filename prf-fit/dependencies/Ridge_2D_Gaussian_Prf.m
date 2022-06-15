function Ridge_2D_Gaussian_Prf(SrfFile, pRFFolder, Roi)
%
% Fits a standard 2D Gaussian pRF model
%   SrfFiles:   Cell array with SamSrf data files (without extension)
%   Roi:        ROI label to restrict analysis 
% Both inputs are optional. If undefined, a dialog is opened for user selection.
%
% This is an example model. Move this file to your parent data folder and
% adapt the model parameters to suit your personal needs and desires.
%
% Rearranged for ERF study, 2019
% samsrf6.05, MATLAB R2016b, SPM12
% ----------------------------------------------------------------------- %
% Standard 2D Gaussian pRF
Model.Prf_Function = @(P,ApWidth) prf_gaussian_rf(P(1), P(2), P(3), ApWidth); % Which pRF model function? 
Model.Name = 'pRF_Gaussian'; % File name to indicate type of pRF model
Model.Param_Names = {'x0'; 'y0'; 'Sigma'}; % Names of parameters to be fitted
Model.Scaled_Param = [1 1 1]; % Which of these parameters are scaled 
Model.Only_Positive = [0 0 1]; % Which parameters must be positive?
Model.Scaling_Factor = 13.5; % Scaling factor of the stimulus space (e.g. eccentricity)
Model.TR = 2; % Repetition time (TR) of pulse sequence
Model.Hrf = []; % HRF file or vector to use (empty = canonical)
Model.Aperture_File = [pRFFolder 'aps_Ridge_mirror']; % Aperture file
% Search grid for coarse fit
Model.Param1 = -1.05 : 0.15 : 1.05; % X0 search grid
Model.Param2 = -1.05 : 0.15 : 1.05; % Y0 search grid
Model.Param3 = 2 .^ (-5.6 : 0.2 : 1); % Sigma search grid
Model.Param4 = 0; % Unused
Model.Param5 = 0; % Unused
% ----------------------------------------------------------------------- %        
% Open dialogs if needed
HomePath = pwd;
% Choose data files
if nargin == 0
    [SrfFile, PathName] = uigetfile('*h_*.mat', 'Choose SamSrf files', 'MultiSelect', 'on');
    if SrfFile ~= 0
        cd(PathName);
    else
        error('No data files selected!');
    end
end
% Choose ROI label
if nargin <= 1
    [Roi, RoiPath] = uigetfile('*.label', 'Choose ROI label');
    if Roi ~= 0 
        Roi = [RoiPath Roi(1:end-6)];
    else
        Roi = '';
    end    
end
% ----------------------------------------------------------------------- %
% Fit pRF model
MapFile = samsrf_fit_prf(Model, SrfFile, pRFFolder, Roi);
% ----------------------------------------------------------------------- %
% Post-processing
load(MapFile);
% Field sign, CMF & smooth
R2_Threshold = 0.05; % R^2 threshold for surface calculations
Eccentricity_Range = [1 Model.Scaling_Factor]; % Eccentricity range for surface calculations
Smoothing_Kernels = [10 3]; % First kernel for field sign & second kernel for everything else
Srf = samsrf_surfcalcs(Srf, Roi, R2_Threshold, Eccentricity_Range, 'S', Smoothing_Kernels, false);
% ----------------------------------------------------------------------- %
% Save again
save(MapFile, 'Srf', 'Model', '-v7.3'); 
% ----------------------------------------------------------------------- %
% Return home
cd(HomePath); 