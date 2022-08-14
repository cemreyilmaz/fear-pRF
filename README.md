# fear-pRF
analysing effect of fearful content on pRF profiles

This repo is used to analyze fMRI data which was collected in the project "The effect of negative emotions on functional activity of human visual cortex"
in Computational & Biological Vision Group, Bilkent University, Turkey. 

The fMRI protocol was based on contracting/shrinking ring and rotating wedge as used in the original pRF paper (Dumoulin&Wandell, 2008). The images from
Nencki Affective Picture Set (NAPS; Riegel et al., 2016) were grouped as neutral and fearful, and the scrambled images was created by block scrambling of 
fearful images. The experiment included three sessions with the pRF stimuli rendered by scrambled, neutral, and fearful images. 

The data was fitted to the pRF model by using SamSrf toolbox on MatLab (Schwarzkopf, de Haas, & Alvarez, 2018). The scritps for the pRF modelling are in
fear-pRF/prf-fit folder. Then, the number of reliable vertices, the pRF sizes and locations were compared by using repeated measures ANOVA follwed by 
multiple comparison test as post-hoc analysis. All the scripts for statistical analyses are in fear-pRF/stats folder. The data was visualized with the 
scripts in fear-pRF/visualisiton.
