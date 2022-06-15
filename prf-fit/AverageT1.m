% AverageT1
% written by Cemre Yilmaz in May 22th, 2019
% MATLAB R2016b, SPM12
% Creates average structural file for three conditions: emo, neutral and sc
% All T1 files must be in FS/subjects directory
% ----------------------------------------------------------------------- %
Subj = 'ERF01'; % subject name
FSFolder = '/Applications/freesurfer/subjects';
mri1 = [FSFolder filesep Subj 'sc_T1.nii']; % Condition1
mri2 = [FSFolder filesep Subj 'neutr_T1.nii']; % Condition2
mri3 = [FSFolder filesep Subj 'emo_T1.nii']; % Condition2
ref = mri1; % Use the first session as reference image
source1 = mri2; % to align the image of Condition2 to the reference
source2 = mri3; % to align the image of Condition2 to the reference
% ----------------------------------------------------------------------- %
% Coregistration: Estimate
% Coregistrations of conditions to scrambled
% ----------------------------------------------------------------------- %
matlabbatch{1}.spm.spatial.coreg.estimate.ref = {ref};
matlabbatch{1}.spm.spatial.coreg.estimate.source = {source1};
spm_jobman('run',matlabbatch); clear matlabbatch

matlabbatch{1}.spm.spatial.coreg.estimate.ref = {ref};
matlabbatch{1}.spm.spatial.coreg.estimate.source = {source2};
spm_jobman('run',matlabbatch); disp('alignments of anatomicals: done'); clear matlabbatch
% ----------------------------------------------------------------------- %
% let's average
% ----------------------------------------------------------------------- %
matlabbatch{1}.spm.util.imcalc.input = {mri1; mri2; mri3};
matlabbatch{1}.spm.util.imcalc.output = Subj;
matlabbatch{1}.spm.util.imcalc.outdir = {FSFolder};
matlabbatch{1}.spm.util.imcalc.expression = '(i1+i2+i3)/3';
matlabbatch{1}.spm.util.imcalc.options.interp = -7;
spm_jobman('run',matlabbatch); disp('averaging: done'); clear matlabbatch