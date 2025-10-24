# ERK activity during early blastocyst development

This repository contains the code necessary to reproduce the simulation and analysis shown in the reasearch article titled: "TBD". 

### Data analysis

All data analysis is done in Python and contained within the folder `data_analysis` in this repository. In there you will find the two files that reproduce the results shown in **Figure 5** and **Supplementary Figure 1**. The files that reproduce these results are: 

- `remove_debris.ipynb`: Jupyter notebook with the instructions to compute the size threshold separating debris from nuclei.
- `quantification.ipynb`: Jupyter notebook with the instructions to segment nuclei and quantify the fluorescence using the resulting masks. 

All segmentations rely on `StarDist` and the `tilesegment` package. See [https://github.com/dsb-lab/tilesegment](https://github.com/dsb-lab/tilesegment) for details about installation.

### Simulations

Modeling is performed using mainly the Julia programming language and contained in the `models/models` folder within this repository. A results folder is predefined in `models/results` so the user can save there the resulting figures of the simulations that reproduce the results of the article. In this work, we study 9 models, each contained in its specific folder in `models/models/MODEL_NAME`. The models shown are the following. 

- `bistability`: Simple bistability circuit shown in **Supplementary Figure 2**. 
- `tristability`: Simple tristability circuit shown in **Supplementary Figure 2**.

These two first models are written in Python. The rest are all written in Julia. 

- `saiz_hadjantonakis_2020`: 3d version of the model shown in [Saiz et al. (2020)]( https://doi.org/10.7554/eLife.56079). The model results are shown in **Figure 2** of this work.
- `saiz_hadjantonakis_2020_FGFko`: *in silico* FGF knock out of the previous model. The model results are shown in **Figure 3** and **Supplementary Figure 3**.
- `ERK1_saiz_extension`: Direct extension of the previous two models which includes explicit ERK dyncamics. The model results are shown in **Figure 4**.
- `ERK2_additive`: Alternative ERK model with additive inhibition of NANOG by GATA6 and ERK. The model results are shown in **Figure 6**.
- `ERK2_multiplicative`: Alternative ERK model with multiplicative inhibition of NANOG by GATA6 and ERK. The model results are shown in **Figure 7**.
- `ERK3_FGFR2`: Extension of the multiplicative model which includes explicit dynamics of the FGFR2. The model results are shown in **Figure 8**.
- `ERK3_homeostasis`: Analysis of signaling range, role of initial conditions and model homeostasis of the previous model. The model results are shown in **Figure 9** and **Figure 10**.