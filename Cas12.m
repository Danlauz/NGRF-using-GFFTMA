clear all
cas=12;

%%% The first variable is a Non-Gaussian Random Field (Positive Rank Asymmetry) %%%  
%%% The second variable is a Non-Gaussian Random Field (Positive Rank Asymmetry by linearly changing the shape parameter) %%%   
%%% The third variable is a Non-Gaussian Random Field (Directional Asymmetry)%%%  
%%% No Cross-Rank Asymmetry %%%  

%%% Model Admissible %%%

%% Simulation parameters
nx=500; ny=500; nbsimul=200; nbfam = 21;
dx= 1; dy=1;
seed = 2515235;

%% Multivariate covariance models (Non-Linear Model of Coregionalization, a.k.a symmetrical covariances)
model{1,1} = [2 120/3]    ; model{1,2}=[7 100/1.8634] ; model{1,3}=[7 70/2.5238];
model{2,1} = model{1,2}  ; model{2,2}=[8 100/4]      ; model{2,3}=[8 70/5.3684];
model{3,1} = model{1,3}  ; model{3,2}=model{2,3}    ; model{3,3}=[8 100/4];

c{1,1} = 1; c{1,2}=0.6;  c{1,3}=0.4;
c{2,1} = c{1,2}; c{2,2}=1; c{2,3}=0.4;
c{3,1} = c{1,3}; c{3,2}=c{2,3}; c{3,3}=1;

nu{1,1} = 1; nu{1,2}=2;  nu{1,3}=1.5;
nu{2,1} = nu{1,2}; nu{2,2}=1; nu{2,3}=2;
nu{3,1} = nu{1,3}; nu{3,2}=nu{2,3}; nu{3,3}=1;

%% Parameters of the familly
param.rangeX = cell(3, 3); param.rangeX{1,1} = [10 230]/3;   
param.nu = cell(3, 3); param.nu{2,2} = [0.5 1.5];
param.c = cell(3, 3);
param.dir = cell(3, 1); param.dir{1}=[0 0 ; 0 0]; param.dir{2}=[0 0 ; 0 0]; param.dir{3}=[0 0 ; 0 40];
param.Gaus = [false false false];
param.Below = [true true true];
param.HD = [];
%% Model admissibility
s=[0.0005:0.0001:0.01 0.02:0.001:0.1 0.2:0.01:7 7.2:0.1:100 101:1:1000]'; % set of frenquencies to check
icode = NG_tasc3d(model,c,nu,s, param, nbfam);
disp(['Family is admissible ', num2str(mean(icode) * 100), '% of the time.']);
%% Multivariate Non-Gaussian Random Fields 
tic
datasimNG = GFFTMA_NG(model, c, nu, param, seed, nbsimul, nbfam, nx, dx, ny, dy);
tNG = toc;

%% Figures
Figures