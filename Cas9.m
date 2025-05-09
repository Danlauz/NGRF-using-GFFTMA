clear all
cas=9;

%%% The first variable is a Non-Gaussian Random Field (Positive Rank Asymmetry with linear change of the anisotropy) %%%  
%%% The second variable is a Gaussian Random Field %%%  
%%% Cross-Rank Asymmetry between first and second variable %%%  

%%% Model Admissible %%%

%% Simulation parameters
nx=500; ny=500; nbsimul=200; nbfam = 21;
dx= 1; dy=1;
seed = 2515235;

%% Multivariate covariance models (Non-Linear Model of Coregionalization, a.k.a symmetrical covariances)
model{1,1} = [2 120/3 120/3 120/3 0 0 0]; model{1,2}=[7 100/1.8634 100/1.8634 100/1.8634 0 0 0];
model{2,1} = model{1,2}  ; model{2,2}=[8 100/4 100/4 100/4 0 0 0];
c{1,1} = 1; c{1,2}=0.6; c{2,1} = c{1,2}  ; c{2,2}=1;
nu{1,1} = 1; nu{1,2}=2; nu{2,1} = nu{1,2}  ; nu{2,2}=1;

%% Parameters of the familly
param.rangeX = cell(2, 2); param.rangeX{1,1} = [30 210]/3; 
param.rangeY = cell(2, 2);
param.rotX = cell(2, 2); param.rotX{1,1} = [0 180]; 
param.nu = cell(2, 2);
param.c = cell(2, 2);
param.dir = cell(2, 2); 
param.Gaus = [false true];
param.Below = [1 1];
param.HD = [];
%% Model admissibility
s=[]'; % set of frenquencies to check
icode = NG_tasc3d(model,c,nu,s, param, nbfam);
disp(['Family is admissible ', num2str(mean(icode) * 100), '% of the time.']);

%% Multivariate Non-Gaussian Random Fields 
tic
datasimNG = GFFTMA_NG(model, c, nu, param, seed, nbsimul, nbfam, nx, dx, ny, dy);
tNG = toc;

%% Figures
Figures