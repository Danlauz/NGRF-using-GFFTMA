clear all
cas=13;
%%% Model Admissible %%%
nvar = 3;
%% Simulation parameters
nx=500; ny=500; nbsimul=100; nbfam = 21;
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
param.dir = cell(3, 1); param.dir{1}=[0 0 ; 0 0]; param.dir{2}=[0 0 ; 0 0]; param.dir{3}=[0 0 ; -20 20];
param.Gaus = [false false false];
param.Below = [true true true];
param.HD = [];
%% Model admissibility
s=[0.0005:0.0001:0.01 0.02:0.001:0.1 0.2:0.01:7 7.2:0.1:100 101:1:1000]'; % set of frenquencies to check
icode = NG_tasc3d(model,c,nu,s, param, nbfam);
disp(['Family is admissible ', num2str(mean(icode) * 100), '% of the time.']);

%% Reference for conditioning
tic
dataRefNG = GFFTMA_NG(model, c, nu, param, seed, 1, nbfam, nx, dx, ny, dy);
tNG = toc/nbsimul;

%% Unconditional Multivariate Non-Gaussian Random Fields
param.HD = [];
tic
datasimNG_UC = GFFTMA_NG(model, c, nu, param, seed*2+2, nbsimul, nbfam, nx, dx, ny, dy);
tNG_UC = toc/nbsimul;

%% Sampling Data Type A 
type = 'A';

for i=1:nvar
    p = haltonset(2,'Skip',1e3,'Leap',1e2);
    nbsample = 50;
    X0 = unique(ceil(net(p,nbsample) *diag([nx ny])),'row');
    LocData = (X0(:,2)-1)*nx + X0(:,1);
    clear p
    param.HD{i} = [X0 dataRefNG(LocData,1,i)];
    %param.HD{i} = [[20 20] 0.42];
end
paramA=param;
% Conditional Multivariate Non-Gaussian Random Fields 
tic
[datasimNG_C_A, ErrA] = GFFTMA_NG(model, c, nu, param, seed*2+2, nbsimul, nbfam, nx, dx, ny, dy);
tNG_C_A = toc/nbsimul;

% Sampling Data Type B 
type = 'B';

for i=1:nvar
    p = haltonset(2,'Skip',1e3,'Leap',1e2);
    nbsample = 250;
    X0 = unique(ceil(net(p,nbsample) *diag([nx ny])),'row');
    LocData = (X0(:,2)-1)*nx + X0(:,1);
    clear p
    param.HD{i} = [X0 dataRefNG(LocData,1,i)];
end
paramB=param;
tic
[datasimNG_C_B, ErrB] = GFFTMA_NG(model, c, nu, param, seed*2+2, nbsimul, nbfam, nx, dx, ny, dy);
tNG_C_B = toc/nbsimul;


% Sampling Data Type C 
type = 'C';

for i=1:nvar
    p = haltonset(2,'Skip',451*i,'Leap',78*i);
    nbsample = 50;
    X0 = unique(ceil(net(p,nbsample) *diag([nx ny])),'row');
    LocData = (X0(:,2)-1)*nx + X0(:,1);
    clear p
    param.HD{i} = [X0 dataRefNG(LocData,1,i)];
end
paramC=param;
tic
[datasimNG_C_C, ErrC] = GFFTMA_NG(model, c, nu, param, seed*2+2, nbsimul, nbfam, nx, dx, ny, dy);
tNG_C_C = toc/nbsimul;


% Sampling Data Type D 
type = 'D';
a = [25 100 250];
for i=1:nvar
    p = haltonset(2,'Skip',451*i,'Leap',78*i);
    nbsample = a(i);
    X0 = unique(ceil(net(p,nbsample) *diag([nx ny])),'row');
    LocData = (X0(:,2)-1)*nx + X0(:,1);
    clear p
    param.HD{i} = [X0 dataRefNG(LocData,1,i)];
end
paramD=param;
tic
[datasimNG_C_D, ErrD] = GFFTMA_NG(model, c, nu, param, seed*2+2, nbsimul, nbfam, nx, dx, ny, dy);
tNG_C_D = toc/nbsimul;

%% Figures
Figures_Cond