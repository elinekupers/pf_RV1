# pf_RV1
Code and software to calculate the contribution of RGC and V1 cortical magnification factor (CMF) to perceptual performance fields effects.

This is a MATLAB code repository (tested on version R2018) and is described in the following 2 publications:

1. **Cortical magnification in human visual cortex parallels task performance around the visual field**\
   Noah C Benson, Eline R Kupers, Antoine Barbot, Marisa Carrasco, Jonathan Winawer\
   (2021). _eLife_ 10:e67685. https://doi.org/10.7554/eLife.67685

2. **Asymmetries around the visual field: From retina to cortex to behavior**\
   Eline R Kupers, Noah C Benson,  Marisa Carrasco, Jonathan Winawer\
   (2022). _PLOS Computational Biology_ 18(1): e1009771. https://doi.org/10.1371/journal.pcbi.1009771

Stored simulations for the 2022 PLOS CB paper can be found on the [Open Science Framework](https://osf.io/ywu5v/)

## Toolbox dependencies
**Github Repositories**
* [isetbio](https://github.com/isetbio/isetbio.git)
* [isetcam](https://github.com/ISET/isetcam.git)
* [JWLOrientedGabor](https://github.com/isetbio/JWLOrientedGabor.git)
* [rgcDisplacementMap](https://github.com/gkaguirrelab/rgcDisplacementMap.git)

**MATLAB add-on dependencies**
* Statistics and machine learning toolbox (for fitcsvm)
* SPHERE3D (for rgcDisplacementMap)

## Folder structure
* scripts                   :   Folder with different scripts and functions to recreate individual paper figures.
* external                  :   External functions from other code repositories.
  
* s_runRGCmodelExamples.m   :   Script with examples of how to run different types of observer model simulations
* pfRV1rootPath             :   Function to set rootpath to folder base
* makeAllFigures.m          :   Script to recreate all data figures from the Benson, Kupers, Barbot, Carrasco, Winawer 2021 paper

\
\
GNU General Public License v3.0 
© ER Kupers 2022
