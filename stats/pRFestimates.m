% Save pRF estimates in reliable vertices for each subjects
% Written by C. Yilmaz, Boyaci Lab, 2019
% MATLAB R2016b, SamSrf6.05
% We use unsmoothed data here
% ----------------------------------------------------------------------- %
clear all;
FSFolder = '/Users/cemreyilmaz/Documents/fear-pRF/data'; % where are subject folders
SubjectList = {'erf01', 'erf02','erf03','erf04','erf05','erf06'}; % List of subject names
Conditions = {'scrambled', 'neutral', 'emotional'}; % list of condiiton names
Rois = {'occ', 'V1', 'V2', 'V2v', 'V2d', 'V3', 'V3v', 'V3d', 'V4', 'V3A'}; % first is always 'occ'
Val = 'Sigma'; % Sigma: pRF size, Beta: BOLD amplitude
Bins = 0:0.75:6.75;
Thrs = [0.05 0.5 6.75]; % R^2 threshold & eccentricity range
EstimateFolder = '/Users/cemreyilmaz/Documents/fear-pRF/data/estimates'; % to store the estimates
addpath(genpath('/Users/cemreyilmaz/Documents/fear-pRF/codes')); % add all the relevant codes to the path
% ----------------------------------------------------------------------- %
T_val = table();
T_center = table();
for i = 1:length(Rois)
    % Empty arrays to store data
    Scrambled = [];
    Neutral = [];
    Emotional = [];
    X_sc = [];
    X_neutral = [];
    X_emo = [];
    Y_sc = [];
    Y_neutral = [];
    Y_emo = [];
    Subjects = {};
    SubjCenter = {};
    Hem_Center = {};
    Hemispheres = {};
    Roi = Rois{i}; % define current ROI
    for sNo = 1:length(SubjectList)
        Subject = SubjectList{sNo}; % define current subject
        SubjFolder = [FSFolder filesep Subject]; % the folder of subject includes the FS recon-all outputs and pRF folders for each condition
        % load pRF estimations
        % Conditions{1}: scrambled
        cd([SubjFolder filesep Conditions{1}]) 
        load('lh_pRF_Gaussian.mat');
        lh_Srf1 = samsrf_expand_srf(Srf);
        lh_rv1 = samsrf_loadlabel([SubjFolder filesep 'ROIs' filesep 'lh_' Roi]);
        lh_Srf1.Values(7:8) = [];
        load('rh_pRF_Gaussian.mat');
        rh_Srf1 = samsrf_expand_srf(Srf);
        rh_rv1 = samsrf_loadlabel([SubjFolder filesep 'ROIs' filesep 'rh_' Roi]);
        rh_Srf1.Values(7:8) = [];
        % Conditions{2}: neutral
        cd([SubjFolder filesep Conditions{2}]) 
        load('lh_pRF_Gaussian.mat');
        lh_Srf2 = samsrf_expand_srf(Srf);
        lh_rv2 = samsrf_loadlabel([SubjFolder filesep 'ROIs' filesep 'lh_' Roi]);
        lh_Srf2.Values(7:8) = [];
        load('rh_pRF_Gaussian.mat');
        rh_Srf2 = samsrf_expand_srf(Srf);
        rh_rv2 = samsrf_loadlabel([SubjFolder filesep 'ROIs' filesep 'rh_' Roi]);
        rh_Srf2.Values(7:8) = [];
        % Conditions{3}: emotional
        cd([SubjFolder filesep Conditions{3}]) 
        load('lh_pRF_Gaussian.mat');
        lh_Srf3 = samsrf_expand_srf(Srf);
        lh_rv3 = samsrf_loadlabel([SubjFolder filesep 'ROIs' filesep 'lh_' Roi]);
        lh_Srf3.Values(7:8) = [];
        load('rh_pRF_Gaussian.mat');
        rh_Srf3 = samsrf_expand_srf(Srf);
        rh_rv3 = samsrf_loadlabel([SubjFolder filesep 'ROIs' filesep 'rh_' Roi]);
        rh_Srf3.Values(7:8) = [];
        % we will use raw data, not the smoothed data
        cd(SubjFolder)
        lh_Srf1.Data = lh_Srf1.Raw_Data; %sc
        rh_Srf1.Data = rh_Srf1.Raw_Data; %sc
        lh_Srf2.Data = lh_Srf2.Raw_Data; %neutral
        rh_Srf2.Data = rh_Srf2.Raw_Data; %neutral
        lh_Srf3.Data = lh_Srf3.Raw_Data; %emo
        rh_Srf3.Data = rh_Srf3.Raw_Data; %emo
        % extract the data from current ROI
        lh_Srf1.Data = lh_Srf1.Data(:,lh_rv1); %sc
        rh_Srf1.Data = rh_Srf1.Data(:,rh_rv1); %sc
        lh_Srf2.Data = lh_Srf2.Data(:,lh_rv2); %neutral
        rh_Srf2.Data = rh_Srf2.Data(:,rh_rv2); %neutral
        lh_Srf3.Data = lh_Srf3.Data(:,lh_rv3); %emo
        rh_Srf3.Data = rh_Srf3.Data(:,rh_rv3); %emo
        % calculate eccentricity
        lh_E1 = sqrt(lh_Srf1.Data(2,:).^2 + lh_Srf1.Data(3,:).^2); %sc
        rh_E1 = sqrt(rh_Srf1.Data(2,:).^2 + rh_Srf1.Data(3,:).^2); %sc
        lh_E2 = sqrt(lh_Srf2.Data(2,:).^2 + lh_Srf2.Data(3,:).^2); %neutral
        rh_E2 = sqrt(rh_Srf2.Data(2,:).^2 + rh_Srf2.Data(3,:).^2); %neutral
        lh_E3 = sqrt(lh_Srf3.Data(2,:).^2 + lh_Srf3.Data(3,:).^2); %emo
        rh_E3 = sqrt(rh_Srf3.Data(2,:).^2 + rh_Srf3.Data(3,:).^2); %emo
        % calculate polar
        lh_P1 = atan2(lh_Srf1.Data(3,:), lh_Srf1.Data(2,:)) / pi * 180; %sc 
        rh_P1 = atan2(rh_Srf1.Data(3,:), rh_Srf1.Data(2,:)) / pi * 180; %sc
        lh_P2 = atan2(lh_Srf2.Data(3,:), lh_Srf2.Data(2,:)) / pi * 180; %neutral
        rh_P2 = atan2(rh_Srf2.Data(3,:), rh_Srf2.Data(2,:)) / pi * 180; %neutral
        lh_P3 = atan2(lh_Srf3.Data(3,:), lh_Srf3.Data(2,:)) / pi * 180; %emo
        rh_P3 = atan2(rh_Srf3.Data(3,:), rh_Srf3.Data(2,:)) / pi * 180; %emo
        % eliminate the data according to the thresholds
        lh_ge = lh_Srf1.Data(1,:) >= Thrs(1) ...
            & lh_Srf2.Data(1,:) >= Thrs(1) ...
            & lh_Srf3.Data(1,:) >= Thrs(1) ...
            & lh_E1 > Thrs(2) & lh_E2 > Thrs(2) & lh_E3 > Thrs(2)  ...
            & lh_E1 < Thrs(3) & lh_E2 < Thrs(3) & lh_E3 < Thrs(3) ...
            & lh_Srf1.Data(4,:) > 0 & lh_Srf2.Data(4,:) > 0 & lh_Srf3.Data(4,:) > 0;
        rh_ge = rh_Srf1.Data(1,:) >= Thrs(1) ...
            & rh_Srf2.Data(1,:) >= Thrs(1) ...
            & rh_Srf3.Data(1,:) >= Thrs(1) ...
            & rh_E1 > Thrs(2) & rh_E2 > Thrs(2) & rh_E3 > Thrs(2) ...
            & rh_E1 < Thrs(3) & rh_E2 < Thrs(3) & rh_E3 < Thrs(3) ...
            & rh_Srf1.Data(4,:) > 0 & rh_Srf2.Data(4,:) > 0 & rh_Srf3.Data(4,:) > 0;
        % Extract the non-zero data in current ROI from Srf1
        switch lower(Val)
            case 'polar'
                lh_Scrambled = lh_P1;
                rh_Scrambled = rh_P1;
            case 'eccentricity'
                lh_Scrambled = lh_E1;
                rh_Scrambled = rh_E1;
            otherwise
                % Left hemisphere
                lh_vs = strfind(lh_Srf1.Values, Val);
                for i = 1:length(lh_Srf1.Values)
                    if isempty(lh_vs{i})
                        lh_vs{i} = 0;
                    end
                end
                lh_xs = strfind(lh_Srf1.Values, 'x0');
                for i = 1:length(lh_Srf1.Values)
                    if isempty(lh_xs{i})
                        lh_xs{i} = 0;
                    end
                end
                lh_ys = strfind(lh_Srf1.Values, 'y0');
                for i = 1:length(lh_Srf1.Values)
                    if isempty(lh_ys{i})
                        lh_ys{i} = 0;
                    end
                end
                lh_vs = cell2mat(lh_vs);
                lh_vs = find(lh_vs);
                lh_xs = cell2mat(lh_xs);
                lh_xs = find(lh_xs);
                lh_ys = cell2mat(lh_ys);
                lh_ys = find(lh_ys);
                lh_Scrambled = lh_Srf1.Data(lh_vs,lh_ge);
                lh_X_sc = lh_Srf1.Data(lh_xs,lh_ge);
                lh_Y_sc = lh_Srf1.Data(lh_ys,lh_ge);
                % Right hemisphere
                rh_vs = strfind(rh_Srf1.Values, Val);
                for i = 1:length(rh_Srf1.Values)
                    if isempty(rh_vs{i})
                        rh_vs{i} = 0;
                    end
                end
                rh_xs = strfind(rh_Srf1.Values, 'x0');
                for i = 1:length(rh_Srf1.Values)
                    if isempty(rh_xs{i})
                        rh_xs{i} = 0;
                    end
                end
                rh_ys = strfind(rh_Srf1.Values, 'y0');
                for i = 1:length(rh_Srf1.Values)
                    if isempty(rh_ys{i})
                        rh_ys{i} = 0;
                    end
                end
                rh_vs = cell2mat(rh_vs);
                rh_vs = find(rh_vs);
                rh_xs = cell2mat(rh_xs);
                rh_xs = find(rh_xs);
                rh_ys = cell2mat(rh_ys);
                rh_ys = find(rh_ys);
                rh_Scrambled = rh_Srf1.Data(rh_vs,rh_ge);
                rh_X_sc = rh_Srf1.Data(rh_xs,rh_ge);
                rh_Y_sc = rh_Srf1.Data(rh_ys,rh_ge);
        end
        % Extract the relevant data from Srf2
        switch lower(Val)
            case 'polar'
                lh_Neutral = lh_P2;
                rh_Neutral = rh_P2;
            case 'eccentricity'
                lh_Neutral = lh_E2;
                rh_Neutral = rh_E2;
            otherwise
                % Left hemisphere
                lh_vs = strfind(lh_Srf2.Values, Val);
                for i = 1:length(lh_Srf2.Values)
                    if isempty(lh_vs{i})
                        lh_vs{i} = 0;
                    end
                end
                lh_xs = strfind(lh_Srf2.Values, 'x0');
                for i = 1:length(lh_Srf2.Values)
                    if isempty(lh_xs{i})
                        lh_xs{i} = 0;
                    end
                end
                lh_ys = strfind(lh_Srf2.Values, 'y0');
                for i = 1:length(lh_Srf2.Values)
                    if isempty(lh_ys{i})
                        lh_ys{i} = 0;
                    end
                end
                lh_vs = cell2mat(lh_vs);
                lh_vs = find(lh_vs);
                lh_xs = cell2mat(lh_xs);
                lh_xs = find(lh_xs);
                lh_ys = cell2mat(lh_ys);
                lh_ys = find(lh_ys);
                lh_Neutral = lh_Srf2.Data(lh_vs,lh_ge);
                lh_X_neutral = lh_Srf2.Data(lh_xs,lh_ge);
                lh_Y_neutral = lh_Srf2.Data(lh_ys,lh_ge);
                % Right hemisphere
                rh_vs = strfind(rh_Srf2.Values, Val);
                for i = 1:length(rh_Srf2.Values)
                    if isempty(rh_vs{i})
                        rh_vs{i} = 0;
                    end
                end
                rh_xs = strfind(rh_Srf2.Values, 'x0');
                for i = 1:length(rh_Srf2.Values)
                    if isempty(rh_xs{i})
                        rh_xs{i} = 0;
                    end
                end
                rh_ys = strfind(rh_Srf2.Values, 'y0');
                for i = 1:length(rh_Srf2.Values)
                    if isempty(rh_ys{i})
                        rh_ys{i} = 0;
                    end
                end
                rh_vs = cell2mat(rh_vs);
                rh_vs = find(rh_vs);
                rh_xs = cell2mat(rh_xs);
                rh_xs = find(rh_xs);
                rh_ys = cell2mat(rh_ys);
                rh_ys = find(rh_ys);
                rh_Neutral = rh_Srf2.Data(rh_vs,rh_ge);
                rh_X_neutral = rh_Srf2.Data(rh_xs,rh_ge);
                rh_Y_neutral = rh_Srf2.Data(rh_ys,rh_ge);
        end
        % Extract the relevant data from Srf3
        switch lower(Val)
            case 'polar'
                lh_Emotional = lh_P3;
                rh_Emotional = rh_P3;
            case 'eccentricity'
                lh_Emotional = lh_E3;
                rh_Emotional = rh_E3;
            otherwise
                % Left hemisphere
                lh_vs = strfind(lh_Srf3.Values, Val);
                for i = 1:length(lh_Srf3.Values)
                    if isempty(lh_vs{i})
                        lh_vs{i} = 0;
                    end
                end
                lh_xs = strfind(lh_Srf3.Values, 'x0');
                for i = 1:length(lh_Srf3.Values)
                    if isempty(lh_xs{i})
                        lh_xs{i} = 0;
                    end
                end
                lh_ys = strfind(lh_Srf3.Values, 'y0');
                for i = 1:length(lh_Srf3.Values)
                    if isempty(lh_ys{i})
                        lh_ys{i} = 0;
                    end
                end
                lh_vs = cell2mat(lh_vs);
                lh_vs = find(lh_vs);
                lh_xs = cell2mat(lh_xs);
                lh_xs = find(lh_xs);
                lh_ys = cell2mat(lh_ys);
                lh_ys = find(lh_ys);
                lh_Emotional = lh_Srf3.Data(lh_vs,lh_ge);
                lh_X_emo = lh_Srf3.Data(lh_xs,lh_ge);
                lh_Y_emo = lh_Srf3.Data(lh_ys,lh_ge);
                % Right hemisphere
                rh_vs = strfind(rh_Srf3.Values, Val);
                for i = 1:length(rh_Srf3.Values)
                    if isempty(rh_vs{i})
                        rh_vs{i} = 0;
                    end
                end
                rh_xs = strfind(rh_Srf3.Values, 'x0');
                for i = 1:length(rh_Srf3.Values)
                    if isempty(rh_xs{i})
                        rh_xs{i} = 0;
                    end
                end
                rh_ys = strfind(rh_Srf3.Values, 'y0');
                for i = 1:length(rh_Srf3.Values)
                    if isempty(rh_ys{i})
                        rh_ys{i} = 0;
                    end
                end
                rh_vs = cell2mat(rh_vs);
                rh_vs = find(rh_vs);
                rh_xs = cell2mat(rh_xs);
                rh_xs = find(rh_xs);
                rh_ys = cell2mat(rh_ys);
                rh_ys = find(rh_ys);
                rh_Emotional = rh_Srf3.Data(rh_vs,rh_ge);
                rh_X_emo = rh_Srf3.Data(rh_xs,rh_ge);
                rh_Y_emo = rh_Srf3.Data(rh_ys,rh_ge);
        end
        % Arrange the data for current subject. Store the extracted data to
        % be concatenated after Subjects for-loop
        sEnd = length(Subjects);
        Subjects = cat(1,Subjects, cell((length(lh_Scrambled)+length(rh_Scrambled)),1));
        Subjects(sEnd+1:end) = {num2str(sNo)};
        lh_CenterSize = [length(lh_X_sc) length(lh_X_emo) length(lh_X_neutral) length(lh_Y_neutral) length(lh_Y_sc) length(lh_Y_emo)];
        lh_X_sc = lh_X_sc(1:min(lh_CenterSize));
        lh_X_neutral = lh_X_neutral(1:min(lh_CenterSize));
        lh_X_emo = lh_X_emo(1:min(lh_CenterSize));
        lh_Y_sc = lh_Y_sc(1:min(lh_CenterSize));
        lh_Y_neutral = lh_Y_neutral(1:min(lh_CenterSize));
        lh_Y_emo = lh_Y_emo(1:min(lh_CenterSize));
        rh_CenterSize = [length(rh_X_sc) length(rh_X_emo) length(rh_Y_sc) length(rh_Y_emo)];
        rh_X_sc = rh_X_sc(1:min(rh_CenterSize));
        rh_X_neutral = rh_X_neutral(1:min(rh_CenterSize));
        rh_X_emo = rh_X_emo(1:min(rh_CenterSize));
        rh_Y_sc = rh_Y_sc(1:min(rh_CenterSize));
        rh_Y_neutral = rh_Y_neutral(1:min(rh_CenterSize));
        rh_Y_emo = rh_Y_emo(1:min(rh_CenterSize));
        X_sc = [X_sc; lh_X_sc'; rh_X_sc'];
        X_neutral = [X_neutral; lh_X_neutral'; rh_X_neutral'];
        X_emo = [X_emo; lh_X_emo'; rh_X_emo'];
        Y_sc = [Y_sc; lh_Y_sc'; rh_Y_sc'];
        Y_neutral = [Y_neutral; lh_Y_neutral'; rh_Y_neutral'];
        Y_emo = [Y_emo; lh_Y_emo'; rh_Y_emo'];
        scEnd = length(SubjCenter);
        SubjCenter = cat(1,SubjCenter, cell((length(lh_X_emo)+length(rh_X_emo)),1));
        SubjCenter(scEnd+1:end) = {num2str(sNo)};
        Hemisphere1 = cell(length(lh_Scrambled),1);
        Hemisphere1(1:end) = {'L'};
        Hemisphere2 = cell(length(rh_Scrambled),1);
        Hemisphere2(1:end) = {'R'};
        Hemispheres = cat(1,Hemispheres, Hemisphere1, Hemisphere2);
        Hem_Center1 = cell(length(lh_X_emo),1);
        Hem_Center1(1:end) = {'L'};
        Hem_Center2 = cell(length(rh_X_emo),1);
        Hem_Center2(1:end) = {'R'};
        Hem_Center = cat(1,Hem_Center, Hem_Center1, Hem_Center2);
        Scrambled = [Scrambled; lh_Scrambled'; rh_Scrambled'];
        Neutral = [Neutral; lh_Neutral'; rh_Neutral'];
        Emotional = [Emotional; lh_Emotional'; rh_Emotional'];
    end
    % Save the data for Sc vs Emo in CSV file. The data will be saved for 2 subjects as:
%     Subjects   Hemispheres   Scrambled   Neutral      Emotional
%     1           L             s1          n1            e1
%     1           L             s2          n2            e2
%     1           L             s3          n3            e3
%     1           R             s1          n1            e1
%     1           R             s2          n2            e2
%     1           R             s3          n3            e3
%     2           L             s1          n1            e1
%     2           L             s2          n2            e2
%     2           L             s3          n3            e3
%     2           R             s1          n1            e1
%     2           R             s2          n2            e2
%     2           R             s3          n3            e3
    ROI = repmat({Roi},length(Subjects),1);
    T = table(Subjects, Hemispheres, ROI, Scrambled, Neutral, Emotional);
    T_val = cat(1,T,T_val);
%     cd(EstimateFolder);
%     save([Val '_' Roi '.mat'], 'Scrambled',  'Neutral', 'Emotional');
%     writetable(T, [Val '_' Roi '.csv']);
%     disp(['Data of ' Val ' were saved in: ' EstimateFolder filesep Val '_' Roi '.csv'])
    Scrambled = sqrt(X_sc.^2+Y_sc.^2);
    Neutral = sqrt(X_neutral.^2+Y_neutral.^2);
    Emotional = sqrt(X_emo.^2+Y_emo.^2);
    Sc_Polar = atan2(Y_sc, X_sc) / pi * 180;
    Neutral_Polar = atan2(Y_neutral, X_neutral) / pi * 180;
    Emo_Polar = atan2(Y_emo, X_emo) / pi * 180;
    Subjects = SubjCenter;
    Hemispheres = Hem_Center;
    T = table(Subjects, Hemispheres, ROI, Scrambled, Neutral, Emotional, Sc_Polar, Neutral_Polar, Emo_Polar);
    T_center = cat(1,T_center,T);
%     save(['pRFcenter_' Roi '.mat'], 'Scrambled',  'Neutral', 'Emotional', 'Sc_Polar', 'Neutral_Polar', 'Emo_Polar');
%     writetable(T_center, ['pRFcenter_' Roi '.csv']);
%     disp(['Data of center locations were saved in: ' EstimateFolder filesep 'pRFcenter_' Roi '.csv'])
end 
cd(EstimateFolder);
writetable(T_val, [Val '.csv']);
writetable(T_center, 'pRFcenter.csv');
% ----------------------------------------------------------------------- %
% Vertex Number for Occipital Lobe
% Define the variables
ValDv = 'Sigma';
ValIv = 'Eccentricity';
Bins = 0:0.75:6.75; % minEccentricity:stepSize:maxEccentricity
Threshold = 0.05;
Mode = 'Mean';
% Create empty arrays to store the data
lh_Scrambled = [];
rh_Scrambled = [];
lh_Neutral = [];
rh_Neutral = [];
lh_Emotional = [];
rh_Emotional = [];
% Let's extract the vertices number for each subject
for sNo = 1:length(SubjectList)
    Subject = SubjectList{sNo};
    % Scrambled
    sc_pRFFolder = [FSFolder filesep Subject filesep Conditions{1}];
    cd(sc_pRFFolder)
    Roi = [FSFolder filesep Subject filesep 'ROIs' filesep 'lh_' Rois{1}]; % Rois{1} = 'occ'
    load('lh_pRF_Gaussian')
    lh_s = vertices(Srf, ValDv, ValIv, Bins, Roi, Threshold, Mode);
    Roi = [FSFolder filesep Subject filesep 'ROIs' filesep 'rh_' Rois{1}]; % Rois{1} = 'occ'
    load('rh_pRF_Gaussian')
    rh_s = vertices(Srf, ValDv, ValIv, Bins, Roi, Threshold, Mode);
    lh_Scrambled = [lh_Scrambled; lh_s];
    rh_Scrambled = [rh_Scrambled; rh_s];
    % Neutral
    neutr_pRFFolder = [FSFolder filesep Subject filesep Conditions{2}];
    cd(neutr_pRFFolder)
    Roi = [FSFolder filesep Subject filesep 'ROIs' filesep 'lh_' Rois{1}]; % Rois{1} = 'occ'
    load('lh_pRF_Gaussian')
    lh_n = vertices(Srf, ValDv, ValIv, Bins, Roi, Threshold, Mode);
    Roi = [FSFolder filesep Subject filesep 'ROIs' filesep 'rh_' Rois{1}]; % Rois{1} = 'occ'
    load('rh_pRF_Gaussian')
    rh_n = vertices(Srf, ValDv, ValIv, Bins, Roi, Threshold, Mode);
    lh_Neutral = [lh_Neutral; lh_n];
    rh_Neutral = [rh_Neutral; rh_n];
    % Emotional
    emo_pRFFolder = [FSFolder filesep Subject filesep Conditions{3}];
    cd(emo_pRFFolder)
    Roi = [FSFolder filesep Subject filesep 'ROIs' filesep 'lh_' Rois{1}]; % Rois{1} = 'occ'
    load('lh_pRF_Gaussian')
    lh_e = vertices(Srf, ValDv, ValIv, Bins, Roi, Threshold, Mode);
    Roi = [FSFolder filesep Subject filesep 'ROIs' filesep 'rh_' Rois{1}]; % Rois{1} = 'occ'
    load('rh_pRF_Gaussian')
    rh_e = vertices(Srf, ValDv, ValIv, Bins, Roi, Threshold, Mode);
    lh_Emotional = [lh_Emotional; lh_e];
    rh_Emotional = [rh_Emotional; rh_e];
end
% Arrange the data. The data will be saved for 3 subjects as:
% Subjects   Hemispheres   Scrambled   Neutral   Emotional
% 1           L             s1          n1        e1
% 2           L             s2          n2        e2
% 3           L             s3          n3        e3
% 1           R             s1          n1        e1
% 2           R             s2          n2        e2
% 3           R             s3          n3        e3
Subjects = 1:length(SubjectList);
Hemisphere1 = cell(length(Subjects),1);
Hemisphere1(1:end) = {'L'};
Hemisphere2 = cell(length(Subjects),1);
Hemisphere2(1:end) = {'R'};
Hemispheres = cat(1, Hemisphere1, Hemisphere2);
Subjects = [Subjects'; Subjects'];
Scrambled = [lh_Scrambled; rh_Scrambled];
Neutral = [lh_Neutral; rh_Neutral];
Emotional = [lh_Emotional; rh_Emotional];
T = table(Subjects, Hemispheres, Scrambled, Neutral, Emotional);
% Save the data in CSV file
cd(EstimateFolder);
writetable(T, ['Vertices_' ValIv '_vs_' ValDv '.csv']);
disp(['Data of reliable vertices number were saved in: ' EstimateFolder filesep 'Vertices_' ValIv '_vs_' ValDv '.csv'])