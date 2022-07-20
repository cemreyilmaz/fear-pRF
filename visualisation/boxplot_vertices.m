% boxplot of number of reliable vertices
% written by C.Yilmaz, Boyaci Lab, 2022
% MATLAB R2022a
clear all; close all;
EstimateFolder = '/Users/cemreyilmaz/Documents/fear-pRF/data/estimates';
FigFolder = '/Users/cemreyilmaz/Documents/fear-pRF/figures/';
cd(EstimateFolder)
data = readtable('Vertices_Eccentricity_vs_Sigma.csv');
% left hemis
boxplot(table2array(data(1:6,3:5)))
ylabel("number of reliable vertices","FontSize",15)
xtickangle(45)
set(gca,'xtick',1:3,'xticklabel',...
    data.Properties.VariableNames(3:5), 'FontSize', 15)
saveas(gca,[FigFolder,'left_vertices.png'])
% right hemis
boxplot(table2array(data(7:12,3:5)))
xtickangle(45)
set(gca,'xtick',1:3,'xticklabel',...
    data.Properties.VariableNames(3:5), 'FontSize', 15)
saveas(gca,[FigFolder,'right_vertices.png'])
