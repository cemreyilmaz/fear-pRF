% Functional Preprocessing & Analysis
% two coregistration is performed, second is for normalization
% for second coreg, averaged structural is used
% written in MATLAB R2016b by O. Batuhan Erkat, arranged by Cemre Yilmaz
% using samsrf6.05
% Last edit: May 22th, 2019
% ----------------------------------------------------------------------- %
clear all; close all; clc
addpath(genpath('.../codes')); % add the path including all the analysis codes
% ----------------------------------------------------------------------- %
% Things to edit for paths and sequence parameters
% ----------------------------------------------------------------------- %
Subject = 'erf01sc'; % Subject ID to analyze
condition = 'scrambled';
AvgSubject = 'erf01'; % one of the session to which the current session will be aligned
RootFolder = '.../raw'; % where my raw dicoms are at
FSFolder = '.../output'; % freesurfer path
Runname = 'Ridge'; % Ridge1, 2, 3 the run folders
Merging=1; % do you need to merge data? Yes[1] No[0]
Averaging=1; % do you need to average data? Yes[1] No[0]
Projecting=1; % do you need to project data onto reconstructed surface? Yes[1] No[0]
NoOfRuns = 12; % no of functional runs
NoOfMeas = 185; % how many image files are there for each run
Dummy = 5; % how many empty scans at the beginning of each run
Overrun = 0; % how many empty scans at the end of each run
% ----------------------------------------------------------------------- %
% No need to edit, mostly initializations, but check if correct just in case
% ----------------------------------------------------------------------- %
FSSubject = [FSFolder filesep Subject]; 
AvgFolder = [FSFolder filesep AvgSubject];
RootSubject = [RootFolder filesep Subject];
SPMFolder = [FSSubject filesep 'spm' filesep]; 
mriFolder = [FSSubject filesep 'mri' filesep];
pRFFolder = [AvgFolder filesep condition filesep]; 
surfFolder = [AvgFolder filesep 'surf' filesep];
CoregFS_T1 = [mriFolder Subject '.nii,1'];
Avg_T1 = [AvgFolder filesep 'mri' filesep AvgSubject '.nii,1'];
Freesurfer_T1 = [AvgFolder filesep 'mri' filesep AvgSubject];
% Empty cell arrays
Runnames = cell(NoOfRuns, 1); 
IMAfolders = cell(NoOfRuns, 1);
Niifolders = cell(NoOfRuns, 1); 
Imagefiles = cell(NoOfMeas, NoOfRuns);
Imagelocs = cell(NoOfMeas, NoOfRuns);
Niifiles = cell(NoOfMeas, NoOfRuns); 
Niilocations = cell(NoOfMeas, NoOfRuns);
afolders = cell(NoOfMeas, NoOfMeas);
bfolders = cell(NoOfMeas,NoOfMeas);
Allpaths = {}; Allfiles = {}; Finallocs = {}; Deletefiles = {};
% ----------------------------------------------------------------------- %
% Loops for initialization of folders and paths
% ----------------------------------------------------------------------- %
if ~exist(SPMFolder); mkdir(SPMFolder); end % create SPMFolder
if ~exist(pRFFolder); mkdir(pRFFolder); end % create pRFFolder
% Folder name array
for i=1:NoOfRuns; Runnames{i,1} = [Runname num2str(i)]; end % list the run names
for j = 1:NoOfRuns % Root subject data folder paths
    IMAfolders{j,1} = [RootSubject filesep Runnames{j,1}];
    Niifolders{j,1} = [SPMFolder Runnames{j,1}];
end
for m = 1:length(Niifolders) % Makes the functional folder names for output (in SPMFolder)
    if ~exist(char(Niifolders(m))); mkdir(char(Niifolders(m))); end
end
for k = 1:NoOfRuns % Takes the image file names and location, cats it to a fullpath for spm
    for s = 1:NoOfMeas
        afolders(s,k) = IMAfolders(k,1); a = IMAfolders(k,1);
        a = char(a); a = dir(a); a = a(3:end,:);
        Imagefiles{s,k} = a(s).name;
        Imagelocs{s,k} = [afolders{s,k} filesep Imagefiles{s,k}];
    end
end
disp('initializations done');
% ----------------------------------------------------------------------- %
% Does conversion job
% ----------------------------------------------------------------------- %
% DICOM import
spm('defaults','fmri'); spm_jobman('initcfg');
for h = 1:NoOfRuns
    matlabbatch{h}.spm.util.import.dicom.data = Imagelocs(:,h); % where dicom files are for each run
    matlabbatch{h}.spm.util.import.dicom.root = 'flat';
    matlabbatch{h}.spm.util.import.dicom.outdir = Niifolders(h); % output dir for each run
    matlabbatch{h}.spm.util.import.dicom.protfilter = '.*';
    matlabbatch{h}.spm.util.import.dicom.convopts.format = 'nii'; % output format
    matlabbatch{h}.spm.util.import.dicom.convopts.icedims = 0;
end
spm_jobman('run',matlabbatch);
clear matlabbatch;
for w = 1:NoOfRuns % Takes the output image file names and location, cats it to full path
    for t = 1:NoOfMeas % You should substract -Dummy here if you already deleted the dummy volumes by hand
        bfolders(t,w) = Niifolders(w,1);
        b = Niifolders(w,1); b = char(b); b = dir(b); b = b(3:end,:);
        Niifiles{t,w} = b(t).name;
        Niilocations{t,w} = [bfolders{t,w} filesep Niifiles{t,w}];
    end
end
disp('DICOM import done');
% ----------------------------------------------------------------------- %
% Erase dummy and overrun volumes
% ----------------------------------------------------------------------- %
if Overrun
    for y=1:NoOfRuns
        delete(Niilocations{1:Dummy,y});
        delete(Niilocations{(NoOfMeas-Overrun+1):end,y});
    end
else
    for y=1:NoOfRuns; delete(Niilocations{1:Dummy,y}); end
end
Niilocations((NoOfMeas-Overrun+1):end,:) = []; Niilocations(1:Dummy,:) = [];
% Add 1 to the end
for e=1:NoOfMeas-Dummy-Overrun
    for o=1:NoOfRuns; Niilocations{e,o} = [Niilocations{e,o} ',1']; end
end
disp('Dummy files are deleted');
% ----------------------------------------------------------------------- %
% SPM Preprocessing
% Realign and unwarp
% ----------------------------------------------------------------------- %
for p=1:NoOfRuns
    matlabbatch{1}.spm.spatial.realignunwarp.data(p).scans = Niilocations(:,p);
    matlabbatch{1}.spm.spatial.realignunwarp.data(p).pmscan = '';
end
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.quality = 0.9;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.sep = 4;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.fwhm = 5;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.rtm = 0;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.einterp = 2;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.ewrap = [0 0 0];
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.weight = '';
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.jm = 0;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.sot = [];
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 4;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.rem = 1;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.noi = 5;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.rinterp = 4;
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.mask = 1;
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.prefix = 'u';
spm_jobman('run',matlabbatch);
clear matlabbatch
disp('realign and unwarp done');%% Delete the old ones that start with f, mean (m) and r, leaves the u ones.
for w = 1:NoOfRuns % Takes the output image file names and location, cats it to full path
    c = Niifolders(w,1);
    c = char(c); c = dir(c); c = c(3:end,:);
    nooffiles = length({c.name});
    for t = 1:nooffiles
        cfolders(t,w)= Niifolders(w,1); Allfiles{t,w} = c(t).name; 
        Finallocs{t,w} = [cfolders{t,w} filesep Allfiles{t,w}];
        if regexp(Allfiles{t,w}, '^f\w*')
            Deletefiles{end+1,1} = Finallocs{t,w}; Finallocs{t,w} = [];
        elseif regexp(Allfiles{t,w}, '^m\w*')
            Deletefiles{end+1,1} = Finallocs{t,w}; Finallocs{t,w} = [];
        elseif regexp(Allfiles{t,w}, '^r\w*')
            Deletefiles{end+1,1} = Finallocs{t,w}; Finallocs{t,w} = [];
        end
    end
end
Finallocs = Finallocs(:); Finallocs(all(cellfun(@isempty,Finallocs),2), :) = [];
for e=1:length(Finallocs); Finallocs{e} = [Finallocs{e} ',1']; end % Add 1 to the end
delete(Deletefiles{:,1}); % delete them
disp('old files are deleted');
% ----------------------------------------------------------------------- %
% Coregistration: Estimate
% Coregistrations of functionals to its own structural image
% ----------------------------------------------------------------------- %
matlabbatch{1}.spm.spatial.coreg.estimate.ref = {CoregFS_T1};
matlabbatch{1}.spm.spatial.coreg.estimate.source = Finallocs(1);
matlabbatch{1}.spm.spatial.coreg.estimate.other = Finallocs(2:end);
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
spm_jobman('run',matlabbatch); disp('coregistration to structural: estimate done'); clear matlabbatch

% ----------------------------------------------------------------------- %
% Normalisation (Coregistration 2)
% Coregistrations of this session to the first session
% ----------------------------------------------------------------------- %
matlabbatch{1}.spm.spatial.coreg.estimate.ref = {Avg_T1};
matlabbatch{1}.spm.spatial.coreg.estimate.source = {CoregFS_T1};
matlabbatch{1}.spm.spatial.coreg.estimate.other = Finallocs(1:end);
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
spm_jobman('run',matlabbatch); disp('coregistration to first session: estimate & write done'); clear matlabbatch
% ----------------------------------------------------------------------- %
% Preprocessing of BOLD for pRF model fitting
% ----------------------------------------------------------------------- %
NoOfRuns = 12; % How many would you like to analyse? If less then
FolderVector={}; %there are only main runs, no hrf or another folder
Merging=1; Averaging=1; Projecting=1; Hemis = {'lh', 'rh'}; %Hemis = {'rh'}; % you can do prf on one hemisphere
spm_jobman('initcfg');
disp('merging niis into run-wise 4D files');
if Merging; MergeData(SPMFolder, pRFFolder, NoOfRuns, FolderVector); end
disp('averaging odd and even runs');
if Averaging; AverageRidge(pRFFolder, NoOfRuns); end
disp('projecting data onto reconstructed surface');
if Projecting; ProjectData(pRFFolder, surfFolder, Freesurfer_T1, Hemis); end
% ----------------------------------------------------------------------- %
% pRF model fitting
% ----------------------------------------------------------------------- %
SrfFullPaths={}; Filenames = dir(char(pRFFolder)); Filenames = Filenames(3:end,:);
for t = 1:length(Filenames)
    if ~isempty(regexp(Filenames(t).name, '^lh\w*')) && ~isempty(regexp(Filenames(t).name, '\w*mat$'))
        SrfFullPaths{end+1,1} = [pRFFolder Filenames(t).name];
   elseif ~isempty(regexp(Filenames(t).name, '^rh\w*')) && ~isempty(regexp(Filenames(t).name, '\w*mat$'))
        SrfFullPaths{end+1,1} = [pRFFolder Filenames(t).name]; 
    end
end
SrfFiles = {};
for a = 1:length(SrfFullPaths)
    [pth SrfFiles{a,1}] = fileparts(char(SrfFullPaths{a,1}));
end
SrfFiles = reshape(SrfFiles, 2, 2); % As there are two hemis, we can do it like this to feed by columns
StartTime=tic; disp('fitting pRF model...');
Rois = {'lh_occ' 'rh_occ'}; % Hemis = {'lh', 'rh'}; 
pRF_DoSubject_Ridge(pRFFolder, FSFolder, surfFolder, SrfFiles, Rois);