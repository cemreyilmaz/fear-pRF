% z-sore functional data
% MATLAB R2016b
% ----------------------------------------------------------------------- %
function result = zscore_mri(data)
% ----------------------------------------------------------------------- %
% zscore the data voxel-by-voxel
for j = 1:size(data,1)
    for l = 1:size(data,2)
        for m = 1:size(data,3)
            Vox = data(j,l,m,:);
            Vox = Vox(1,:);
            Vox = zscore(Vox);
            data(j,l,m,:) = Vox;
        end
    end
end
% ----------------------------------------------------------------------- %
result = data;
end