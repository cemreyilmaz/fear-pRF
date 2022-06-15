% detrend functional data
% samsrf6.05, MATLAB R2016b, SPM12
function result = detrend_mri(data)
% ----------------------------------------------------------------------- %
% detrend the data voxel-by-voxel
for j = 1:size(data,1)
    for l = 1:size(data,2)
        for m = 1:size(data,3)
            Vox = data(j,l,m,:);
            Vox = Vox(1,:);
            Vox = detrend(Vox);
            data(j,l,m,:) = Vox;
        end
    end
end
% ----------------------------------------------------------------------- %
result = data;
end