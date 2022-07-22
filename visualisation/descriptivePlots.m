% Descriptive plots for each ROI - group results
% Before saving the figures, rearrange the location of legend and the
% ylim([]);
% written by C.Yilmaz, Boyaci Lab, 2019
% MATLB R2016b
% ----------------------------------------------------------------------- %
clear all; close all;
EstimateFolder = '/Users/cemreyilmaz/Documents/fear-pRF/data/estimates/';
StatFolder = '/Users/cemreyilmaz/Documents/fear-pRF/results/';
FigFolder = '/Users/cemreyilmaz/Documents/fear-pRF/figures/';
Rois = {'occ', 'V1', 'V2v', 'V2d', 'V3v', 'V3d', 'V4', 'V3A'};
Variables = {'pRFcenter', 'Sigma', 'Beta'};
% ----------------------------------------------------------------------- %
% calculations
% ----------------------------------------------------------------------- %
for roi = 1:length(Rois)
    Roi = Rois{roi};
    % pRF shift
    data = readtable([EstimateFolder Variables{1} '.csv']);
    if ~strcmp(Roi,'occ')
        data = data(strcmp(data.ROI,Roi),:);
    end
    s_shift = data.Scrambled;
    n_shift = data.Neutral;
    e_shift = data.Emotional;
    N = length(s_shift);
    CI95(roi,:) = tinv([0.025 0.975], N-1);
    ScMean_shift(roi,:) = mean(s_shift);
    NeutrMean_shift(roi,:) = mean(n_shift);
    EmoMean_shift(roi,:) = mean(e_shift);
    ScCI95_shift(roi,:) = bsxfun(@times, std(s_shift)/sqrt(N), CI95(roi,:));
    NeutrCI95_shift(roi,:) = bsxfun(@times, std(n_shift)/sqrt(N), CI95(roi,:));
    EmoCI95_shift(roi,:) = bsxfun(@times, std(e_shift)/sqrt(N), CI95(roi,:));
    % pRF size
    data = readtable([EstimateFolder Variables{2} '.csv']);
    if ~strcmp(Roi,'occ')
        data = data(strcmp(data.ROI,Roi),:);
    end
    s_sigma = data.Scrambled;
    n_sigma = data.Neutral;
    e_sigma = data.Emotional;
    N = length(s_sigma);
    CI95(roi,:) = tinv([0.025 0.975], N-1);
    ScMean_sigma(roi,:) = mean(s_sigma);
    NeutrMean_sigma(roi,:) = mean(n_sigma);
    EmoMean_sigma(roi,:) = mean(e_sigma);
    ScCI95_sigma(roi,:) = bsxfun(@times, std(s_sigma)/sqrt(N), CI95(roi,:));
    NeutrCI95_sigma(roi,:) = bsxfun(@times, std(n_sigma)/sqrt(N), CI95(roi,:));
    EmoCI95_sigma(roi,:) = bsxfun(@times, std(e_sigma)/sqrt(N), CI95(roi,:));
    % BOLD
    data = readtable([EstimateFolder Variables{3} '.csv']);
    if ~strcmp(Roi,'occ')
        data = data(strcmp(data.ROI,Roi),:);
    end
    s_beta = data.Scrambled;
    n_beta = data.Neutral;
    e_beta = data.Emotional;
    N = length(s_beta);
    CI95(roi,:) = tinv([0.025 0.975], N-1);
    ScMean_beta(roi,:) = mean(s_beta);
    NeutrMean_beta(roi,:) = mean(n_beta);
    EmoMean_beta(roi,:) = mean(e_beta);
    ScCI95_beta(roi,:) = bsxfun(@times, std(s_beta)/sqrt(N), CI95(roi,:));
    NeutrCI95_beta(roi,:) = bsxfun(@times, std(n_beta)/sqrt(N), CI95(roi,:));
    EmoCI95_beta(roi,:) = bsxfun(@times, std(e_beta)/sqrt(N), CI95(roi,:));
end
x = 1:length(Rois);
% ----------------------------------------------------------------------- %
% Sc vs Emo
% ----------------------------------------------------------------------- %
figure()
scatter(x,ScMean_shift,'filled','MarkerEdgeColor','black', ...
    'MarkerFaceColor','black','LineWidth',0.75)
hold on
scatter(x,EmoMean_shift,'filled','MarkerEdgeColor','red', ...
    'MarkerFaceColor','red','LineWidth',0.75)
errorbar(x,ScMean_shift,ScCI95_shift(:,1),ScCI95_shift(:,2), ...
    'Color','black','LineStyle','None','LineWidth',2)
errorbar(x,EmoMean_shift,EmoCI95_shift(:,1),EmoCI95_shift(:,2), ...
    'Color','red','LineStyle','None','LineWidth',2)
hold off
xlim([0 11])
xticks(x)
xtickangle(45)
set(gca,'xtick',x,'xticklabel',Rois, 'FontSize', 15)
legend({'Scrambled','Emotional','',''}, 'FontSize', 13)
ylabel('pRF center (in visual degrees)', 'FontSize', 15)
%title('Mean locations of pRF center in each retinotopic region', 'FontSize', 18)
saveas(gca,[FigFolder 'pRFcenter_sc_emo.png'])
% ----------------------------------------------------------------------- %
figure()
scatter(x,ScMean_sigma,'filled','MarkerEdgeColor','black', ...
    'MarkerFaceColor','black','LineWidth',0.75)
hold on
scatter(x,EmoMean_sigma,'filled','MarkerEdgeColor','red', ...
    'MarkerFaceColor','red','LineWidth',0.75)
errorbar(x,ScMean_sigma,ScCI95_sigma(:,1),ScCI95_sigma(:,2), ...
    'Color','black','LineStyle','None','LineWidth',2)
errorbar(x,EmoMean_sigma,EmoCI95_sigma(:,1),EmoCI95_sigma(:,2), ...
    'Color','red','LineStyle','None','LineWidth',2)
hold off
xlim([0 11])
xticks(x)
xtickangle(45)
set(gca,'xtick',x,'xticklabel',Rois, 'FontSize', 15)
legend({'Scrambled','Emotional','',''}, 'FontSize', 13)
ylabel('pRF size (in visual degrees)', 'FontSize', 15)
%title('Mean pRF sizes in each retinotopic region', 'FontSize', 18)
saveas(gca,[FigFolder 'pRFsize_sc_emo.png'])
% ----------------------------------------------------------------------- %
figure()
scatter(x,ScMean_beta,'filled','MarkerEdgeColor','black', ...
    'MarkerFaceColor','black','LineWidth',0.75)
hold on
scatter(x,EmoMean_beta,'filled','MarkerEdgeColor','red', ...
    'MarkerFaceColor','red','LineWidth',0.75)
errorbar(x,ScMean_beta,ScCI95_beta(:,1),ScCI95_beta(:,2), ...
    'Color','black','LineStyle','None','LineWidth',2)
errorbar(x,EmoMean_beta,EmoCI95_beta(:,1),EmoCI95_beta(:,2), ...
    'Color','red','LineStyle','None','LineWidth',2)
hold off
xlim([0 11])
xticks(x)
xtickangle(45)
set(gca,'xtick',x,'xticklabel',Rois, 'FontSize', 15)
legend({'Scrambled','Emotional','',''}, 'FontSize', 13)
ylabel('Amplitude of BOLD (estimated by pRF method)', 'FontSize', 15)
%title('Mean pRF sizes in each retinotopic region', 'FontSize', 18)
saveas(gca,[FigFolder 'BOLD_sc_emo.png'])
% ----------------------------------------------------------------------- %
% Neutr vs Emo
% ----------------------------------------------------------------------- %
figure()
scatter(x,NeutrMean_shift,'filled','MarkerEdgeColor','blue', ...
    'MarkerFaceColor','blue','LineWidth',0.75)
hold on
scatter(x,EmoMean_shift,'filled','MarkerEdgeColor','red', ...
    'MarkerFaceColor','red','LineWidth',0.75)
errorbar(x,NeutrMean_shift,NeutrCI95_shift(:,1),NeutrCI95_shift(:,2), ...
    'Color','blue','LineStyle','None','LineWidth',2)
errorbar(x,EmoMean_shift,EmoCI95_shift(:,1),EmoCI95_shift(:,2), ...
    'Color','red','LineStyle','None','LineWidth',2)
hold off
xlim([0 11])
xticks(x)
xtickangle(45)
set(gca,'xtick',x,'xticklabel',Rois, 'FontSize', 15)
legend({'Neutral','Emotional','',''}, 'FontSize', 13)
%ylabel('pRF center (in visual degrees)', 'FontSize', 15)
%title('Mean locations of pRF center in each retinotopic region', 'FontSize', 18)
saveas(gca,[FigFolder 'pRFcenter_neutr_emo.png'])
% ----------------------------------------------------------------------- %
figure()
scatter(x,NeutrMean_sigma,'filled','MarkerEdgeColor','blue', ...
    'MarkerFaceColor','blue','LineWidth',0.75)
hold on
scatter(x,EmoMean_sigma,'filled','MarkerEdgeColor','red', ...
    'MarkerFaceColor','red','LineWidth',0.75)
errorbar(x,NeutrMean_sigma,NeutrCI95_sigma(:,1),NeutrCI95_sigma(:,2), ...
    'Color','blue','LineStyle','None','LineWidth',2)
errorbar(x,EmoMean_sigma,EmoCI95_sigma(:,1),EmoCI95_sigma(:,2), ...
    'Color','red','LineStyle','None','LineWidth',2)
hold off
xlim([0 11])
xticks(x)
xtickangle(45)
set(gca,'xtick',x,'xticklabel',Rois, 'FontSize', 15)
legend({'Neutral','Emotional','',''}, 'FontSize', 13)
%ylabel('pRF size (in visual degrees)', 'FontSize', 15)
%title('Mean pRF sizes in each retinotopic region', 'FontSize', 18)
saveas(gca,[FigFolder 'pRFsize_neutr_emo.png'])
% ----------------------------------------------------------------------- %
figure()
scatter(x,NeutrMean_beta,'filled','MarkerEdgeColor','blue', ...
    'MarkerFaceColor','black','LineWidth',0.75)
hold on
scatter(x,EmoMean_beta,'filled','MarkerEdgeColor','red', ...
    'MarkerFaceColor','red','LineWidth',0.75)
errorbar(x,NeutrMean_beta,NeutrCI95_beta(:,1),NeutrCI95_beta(:,2), ...
    'Color','blue','LineStyle','None','LineWidth',2)
errorbar(x,EmoMean_beta,EmoCI95_beta(:,1),EmoCI95_beta(:,2), ...
    'Color','red','LineStyle','None','LineWidth',2)
hold off
xlim([0 11])
xticks(x)
xtickangle(45)
set(gca,'xtick',x,'xticklabel',Rois, 'FontSize', 15)
legend({'Neutral','Emotional','',''}, 'FontSize', 13)
%ylabel('Amplitude of BOLD (estimated by pRF method)', 'FontSize', 15)
%title('Mean pRF sizes in each retinotopic region', 'FontSize', 18)
saveas(gca,[FigFolder 'BOLD_neutr_emo.png'])
