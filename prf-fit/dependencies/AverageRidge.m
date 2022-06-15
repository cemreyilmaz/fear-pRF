function AverageRidge(pRFFolder, NoOfRuns)
% Rearranged for ERF study, 2019
% samsrf6.05, MATLAB R2016b, SPM12
BatchNo=1;
if mod(NoOfRuns./2,2)
    error('Odd number of runs - I won''t average! Please check!!')
end
% ----------------------------------------------------------------------- %
% let's do it
OddRunNos=1:2:NoOfRuns;
EvenRunNos=2:2:NoOfRuns;
disp('normalizing odds');
%read in data and normalise
for iRun=1:length(OddRunNos)
    CurrOdd=spm_read_vols(spm_vol([pRFFolder 'Ridge' num2str(OddRunNos(iRun)) '_4D.nii']));
    CurrOdd=detrend_mri(CurrOdd);%detrend time-series
    CurrOdd=zscore_mri(CurrOdd);%z-score time-series
    OddRuns(:,:,:,:,iRun)=CurrOdd;
    disp(['done: ' num2str(iRun) 'odd']);
end
% ----------------------------------------------------------------------- %
disp('reading template header');
%read in template hdr
V=spm_vol([pRFFolder 'Ridge1_4D.nii']);
Hdr=V(1);
Hdr.dt(1)=64;%save as double
OddHdr=Hdr;
EvenHdr=Hdr;
% ----------------------------------------------------------------------- %
disp('averaging odds');
% average across runs
OddAvg = squeeze(mean(OddRuns,5));
PrevFolder=pwd;
cd(pRFFolder);
% ----------------------------------------------------------------------- %
disp('writing back to 3D');
%now write (annoyingly it seems we have to revert to 3D for that and than merge again)
for iVol=1:size(OddAvg,4)
    OddHdr.fname=['OddAverage' num2str(iVol) '.nii'];
    spm_write_vol(OddHdr, squeeze(OddAvg(:,:,:,iVol)));
    end%for iVol
% ----------------------------------------------------------------------- %
disp('merging again');
%now merge stuff again...
for iVol=1:size(OddAvg,4)
    OddEPIs{iVol}=[pwd filesep 'OddAverage' num2str(iVol) '.nii'];
end%for iVol
spm_file_merge(OddEPIs, [pRFFolder 'FirstOddAvg_4D.nii']);
%finally, clean up
delete('OddAverage*.nii');
% ----------------------------------------------------------------------- %
disp('normalizing evens');
%read in data and normalise
for iRun=1:length(EvenRunNos)
    CurrEven=spm_read_vols(spm_vol([pRFFolder 'Ridge' num2str(EvenRunNos(iRun)) '_4D.nii']));
    CurrEven=detrend_mri(CurrEven);
    CurrEven=zscore_mri(CurrEven);
    EvenRuns(:,:,:,:,iRun)=CurrEven;
    disp(['done: ' num2str(iRun) 'even']);
end
% ----------------------------------------------------------------------- %
disp('averaging evens');
%average across runs
EvenAvg = squeeze(mean(EvenRuns,5));
PrevFolder=pwd;
cd(pRFFolder);
% ----------------------------------------------------------------------- %
disp('writing back to 3D');
%now write (annoyingly it seems we have to revert to 3D for that and than merge again)
for iVol=1:size(OddAvg,4)
    EvenHdr.fname=['EvenAverage' num2str(iVol) '.nii'];
    spm_write_vol(EvenHdr, squeeze(EvenAvg(:,:,:,iVol)));
end%for iVol
% ----------------------------------------------------------------------- %
disp('merging again');
%now merge stuff again...
for iVol=1:size(OddAvg,4)
    EvenEPIs{iVol}=[pwd filesep 'EvenAverage' num2str(iVol) '.nii'];
end%for iVol
spm_file_merge(EvenEPIs, [pRFFolder 'ThenEvenAvg_4D.nii']);
%finally, clean up
delete('EvenAverage*.nii');
delete([pRFFolder 'Ridge*.*']);
cd(PrevFolder);
end