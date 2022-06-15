% Structural preprocessing
% written in MATLAB R2016b by O. Batuhan Erkat, arranged by Cemre Yilmaz
% Last edit: May 22th, 2019
% ----------------------------------------------------------------------- %
clear all; close all; clc
addpath(genpath('/Applications/MATLAB/toolbox/SamSrf/analysisTools')) % add the path including all the analysis codes
% ----------------------------------------------------------------------- %
% Things to edit for your PC
% ----------------------------------------------------------------------- %
RootFolder = '/Applications/freesurfer/subjects/raw'; % where my raw dicoms are at
Subject = 'ERF01sc'; % subject initials, folder name of raw and procd data
FSFolder = '/Applications/freesurfer/subjects'; % freesurfer path
T1NoOfMeas = 176; % how many images do you have in T1 scan
% ----------------------------------------------------------------------- %
% No need to edit, mostly initializations, but check if correct just in case
% ----------------------------------------------------------------------- %
% RootSubject = [RootFolder filesep Subject];
T1root = [RootSubject filesep 'T1']; % Raw T1 files
T1out = [RootSubject filesep 'T1out']; % Imported T1 as single nii
Runnames = {}; Outfolders = {};
% ----------------------------------------------------------------------- %
% Initializations for DICOM import only for T1
% ----------------------------------------------------------------------- %
% if ~exist(T1out) % Makes the RootSubject folder for output
    mkdir(T1out)
end
T1outcell = cellstr(T1out);
r = T1root; r = char(r); r = dir(r); r = r(3:end,:); % list all the files in ../Subject/T1/..
T1files = {}; T1paths = {}; T1locs = {};
for n = 1:T1NoOfMeas % rearrange the files as one-column list
    T1files{n} = r(n).name;
    T1paths{n} = T1root;
    T1locs{n,1} = [T1paths{n} filesep T1files{n}];
end
% ----------------------------------------------------------------------- %
% Import T1
% ----------------------------------------------------------------------- %
spm('defaults','fmri'); spm_jobman('initcfg');
matlabbatch{1}.spm.util.import.dicom.data = T1locs; % where dicom files are for each run
matlabbatch{1}.spm.util.import.dicom.root = 'flat';
matlabbatch{1}.spm.util.import.dicom.outdir = T1outcell; % output dir for each run
matlabbatch{1}.spm.util.import.dicom.protfilter = '.*';
matlabbatch{1}.spm.util.import.dicom.convopts.format = 'nii'; % output format
matlabbatch{1}.spm.util.import.dicom.convopts.icedims = 0;
spm_jobman('run',matlabbatch);
clear matlabbatch
disp('DICOM import for T1 is done');
% ----------------------------------------------------------------------- %
% Rename the output T1, and relocate it for Freesurfer recon-all
% ----------------------------------------------------------------------- %
e = T1out; e = char(e); e = dir(e); e = e(3:end,:);
e_in = [T1out filesep e(1).name];
T1out_name = [Subject '_T1.nii'];
e_out = [FSFolder filesep T1out_name];
movefile(e_in, e_out); % Renames and moves the file
rmdir(T1out); % delete the empty directory
disp('T1 file is relocated for recon-all on Freesurfer');
% Now you have raw files renamed in your root folder and T1 at Freesurfer
% folder waiting for you to recon-all.