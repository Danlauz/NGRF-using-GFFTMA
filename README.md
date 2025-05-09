# Non-Gaussian multivariate random fields using Generalized Fast Fourier Transform Moving Average

This repository contains the computer codes associated with the manuscript: *A Generalized FFTMA Approach to Simulate Multivariate Non-Gaussian Random Fields*.

## PROGRAMS
- `Cas1.m to Cas13.m` : Perform the synthetic exemple of the manuscript. Return the image in the subfolder Figures.
- `Cas14.m` : Run the real case study of the smelter.
- `GFFTMA.m` : Adapted code from Liang & al. (2016) to simulation Multivarite Gaussian Random Fields.
- `GFFTMA_NG.m` : Code to simulation Multivarite Non-Gaussian Random Fields.
- `RunningCases` : Automaticaly run the desired folder Cas1.m to Cas14.m

## Example
The folder `Functions` contains all the code related to the results presented in the manuscript, organized by number. Figures from the manuscript are stored in `.png` format in their respective folders, and datasets are also included.

Non-Gaussian random fields and synthetic datasets are generated using the `GFFTMA.m` algorithm provided in this folder. The script `covar.m` is used to generate the covariance matrix, while `grille2.m` and `grille3.m` create 2D and 3D grids, respectively. The function `ECDF.m` computes the empirical cumulative distribution using ranking.

Data set from the smelter are in the excel file or can be retrive in Henderson et al. (2002). 

`ApplyMask.m`, `ll2utm.m`, `utm2ll.m`, `saveGridData.m` are helper function for data transformation or 

`anamor_multi.m` performs Gaussian anamorphosis independently on each variable.
`Figures.m`, `Figures_BaseCase.m`, `Figures_Cond.m` are helper function to reproduce the figures on the manuscript.

`GeoStatFFT.m`. This function computes direct- and cross- spatial statistics in `nD` for up to `nvar` variables.
`GeoStatFFT_ndir.m`. This function post-processes the output of `GeoStatFFT` to compute experimental directional or omnidirectional direct- and cross- spatial statistics.
Available Spatial Statistics:
- **1**: Variograms and cross-variograms
- **2**: Covariograms and cross-covariograms
- **3**: Variograms and pseudo-cross-variograms
- **4**: Centered covariances and cross-covariances (mean computed for the entire field)
- **5**: Non-centered covariances and cross-covariances (bivariate probabilities, for categorical data)
- **6**: Transiograms (for categorical data)
- **7**: Non-ergodic transiograms (for categorical data)
- **8**: Directional asymmetry and cross-directional asymmetry (Bárdossy & Hörning, 2017)
- **9**: Rank asymmetry and cross-rank asymmetry (Guthke, 2013)
- **10**: Rank correlation and cross-rank correlation (Bárdossy & Hörning, 2017)
- **11**: Third-order cumulant of a zero-mean random function (Dimitrakopoulos et al., 2010)

#### Reference:
- Liang, M., Marcotte, D., and Shamsipour, P. (2016). *Simulation of non-linear coregionalization models by FFTMA*. Computers & Geosciences 89, 220–231. doi:10.1016/j.cageo.2016.01.005
- Henderson, P. J., Knight, R. D., and McMartin, I. (2002). *Geochemistry of soils within a 100 km radius of the Horne Cu smelter, Rouyn-Noranda, Québec*. doi:10.4095/213097678
- Lauzon, D. and Hörning, S. (2025) *Efficient Computation on Large Regular Grids of High-Order Spatial Statistics via Fast Fourier Transform*. Computers & Geosciences, 198, 105878. doi:10.1016/j.cageo.2025.105878
- Lauzon, D., Hörning, S, Bardossy, A. *A Generalized FFTMA Approach to Simulate Multivariate Non-Gaussian Random Fields*. (In review)
---
