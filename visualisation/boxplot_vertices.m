% boxplot of number of reliable vertices
% written by C.Yilmaz, Boyaci Lab, 2022
% MATLAB R2022a
clear all; close all;
EstimateFolder = '/Users/cemreyilmaz/Documents/fear-pRF/data/estimates';
FigFolder = '/Users/cemreyilmaz/Documents/fear-pRF/figures/';
cd(EstimateFolder)
data = readtable('Vertices_Eccentricity_vs_Sigma.csv');
% left hemis
boxplot([data.Neutral(1:6),data.Scrambled(1:6),data.Emotional(1:6)])
ylabel("number of reliable vertices","FontSize",15)
xtickangle(45)
set(gca,'xtick',1:3,'xticklabel',...
    {'Neutral','Scrambled','Emotional'}, 'FontSize', 15)
yt = get(gca, 'YTick');
axis([xlim    0  ceil(max(yt)*1.2)])
hold on
plot([2.1 2.9], [1 1]*max(yt)*1.1, '-k',  mean([2.1 2.9]), max(yt)*1.15, '*k')
plot([1.1 1.9], [1 1]*max(yt)*1.1, '-k',  mean([1.1 1.9]), max(yt)*1.15, '*k')
hold off
saveas(gca,[FigFolder,'left_vertices.png'])
% right hemis
boxplot([data.Neutral(7:12),data.Scrambled(7:12),data.Emotional(7:12)])
xtickangle(45)
set(gca,'xtick',1:3,'xticklabel',...
    {'Neutral','Scrambled','Emotional'}, 'FontSize', 15)
yt = get(gca, 'YTick');
axis([xlim    0  ceil(max(yt)*1.2)])
hold on
plot([2.1 2.9], [1 1]*max(yt)*1.1, '-k',  mean([2.1 2.9]), max(yt)*1.15, '*k')
plot([1.1 1.9], [1 1]*max(yt)*1.1, '-k',  mean([1.1 1.9]), max(yt)*1.15, '*k')
hold off
saveas(gca,[FigFolder,'right_vertices.png'])
