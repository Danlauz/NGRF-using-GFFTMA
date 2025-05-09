function [ZBest, ErrBest, UBest] = GFFTMA_NG(model, c, nu, param, seed, nbsimul, nbfam ,nx, dx, ny, dy, nz, dz)

%%%%%%%% WARNING - ANISOTROPIES ONLY WORKS IN 2D %%%%%%%%%
%%%%%%%% Memory storage is not optimized %%%%%%%%%

% This code calls GFFTMA multiple times to generate families of Gaussian
% fields that share similar parameters (i.e., shape, range, anisotropies).
% It is assumed that the isotropic case structure of GFFTMA is always used.
% Afterward, truncation rules are applied to generate non-Gaussian random fields.

% NG_GFFTMA: Function to simulate non-Gaussian non-LMC (but symmetrical) models using FFT-MA
%
% Syntax:
%   datasim = NG_GFFTMA(model, seed, nbsimul, vsiz, nx, dx, ny, dy, nz, dz)
%
% Description:
%   Simulates spatial fields using FFT-MA. The method supports 1D, 2D, or 3D simulations.
%   - 1D: Specify nx and dx.
%   - 2D: Specify nx, dx, ny, and dy.
%   - 3D: Specify all parameters (nx, dx, ny, dy, nz, dz).
%
% Inputs:
%   model   : nvar x nvar cell array of structures defining the covariance models.
%             - Each cell {i, j} contains structures for cross-covariance
%               between variables i and j.
%             - 1D or Isotropic case : [model type, range].
%             - 2D Anisotropic case  : [model type, range1, range2, rotx].
%             - 3D Anisotropic case  : [model type, range1, range2, range3, rotx, roty, rotz].
%   c       : nvar x nvar cell array of sill defining the covariance models.
%   nu      : nvar x nvar cell array of shape parameters defining the covariance models (i.e., for Matern covariances).
%   param   : struc. of parameters to create the family
%             - Each struct contains array of paramters for each structure of familly [minRange, maxRange]
%                   - minRange: Minimum range of the family
%                   - maxRange: Maximum range of the family
%   seed    : Random seed for reproducibility.
%   nbsimul : Number of realizations to generate.
%   nbfam   : Number of families to generate
%   nx, dx  : Number of grid points and step size along x-axis.
%   ny, dy  : (Optional) Number of grid points and step size along y-axis.
%   nz, dz  : (Optional) Number of grid points and step size along z-axis.
%
% Outputs:
%   datasim : Cell array of size [nbsimul, nvar]. Each cell contains a vector of
%             simulated values for each variable in the grid.
%
% Author:
%   D. Lauzon, Mai 2025

% Determine simulation dimensionality
if nargin == 9
    cas = 1; % 1D simulation
    ny = 1; nz = 1; dy = 1; dz = 1;
elseif nargin == 11
    cas = 2; % 2D simulation
    nz = 1; dz = 1;
elseif nargin == 13
    cas = 3; % 3D simulation
else
    error('Check input parameters: Incorrect number of inputs.');
end

% Initialize random seed
rng('default');
rng(round(seed));

% Number of variables
nvar = size(model, 1);

% Enlarge grid size to avoir aliasing
Nx=nx*2+1; Ny=ny*2+1; Nz=nz*2+1;

% Compute grid midpoints
Nx2 = floor(Nx / 2);
Ny2 = floor(Ny / 2);
Nz2 = floor(Nz / 2);

% Generate grid coordinates
if cas == 1
    x0 = (-(Nx2 + 1) * dx : dx : Nx2 * dx)';
    x0r = 0;
elseif cas == 2
    x0 = grille2(-(Nx2 + 1) * dx, Nx2 * dx, dx, -(Ny2 + 1) * dy, Ny2 * dy, dy);
    x0r = [0 0];
elseif cas == 3
    x0 = grille3(-(Nx2 + 1) * dx, Nx2 * dx, dx, -(Ny2 + 1) * dy, Ny2 * dy, dy, -(Nz2 + 1) * dz, Nz2 * dz, dz);
    x0r = [0 0 0];
end

% Parameters tau that ensures a continuous transition of C(h,tau)
if nbfam ==1
    tau = 0.5;
else
    tau = 0 : 1/(nbfam-1) : 1;
end

% Maximum number of iterations for the Monte Carlo conditioning step (not perfectly program).
ittmax = 0;

% Initialization of the NonGaussian random fields
ZBest = zeros(nx*ny*nz,nbsimul,nvar);
ErrBest = nan(ittmax, nbsimul);

% Compute the spectral representation of the N-LCM (we store all required frequency matrices).
% This is not the most memory-efficient approach, but it reduces computation time as we compute them only once.
G = cell(nbfam,1);
% Initialize directional asymmetry
vx = zeros(nbfam,nvar); vy = zeros(nbfam,nvar); vz = zeros(nbfam,nvar);
for nf = 1 : nbfam
    %Local variable for parfor loop (if used)
    model_local=model;     nu_local = nu;     c_local=c;
    % Initialisation of temporary variable
    model_temp=cell(nvar, nvar); c_temp=cell(nvar, nvar); nu_temp=cell(nvar, nvar);
    for i = 1:nvar
        for j = i:nvar
            % Update covariance model parameters using a helper function
            % (this computes the covariance model for family number nf)
            model_temp{i, j} = updateParameters(model_local{i, j}, param, tau(nf), i, j, nbfam, {'rangeX', 2; 'rangeY', 3; 'rotX', 5},nf);
            nu_temp{i, j} = updateParameters(nu_local{i, j}, param, tau(nf), i, j, nbfam, {'nu', []},nf);
            c_temp{i, j} = updateParameters(c_local{i, j}, param, tau(nf), i, j, nbfam, {'c', []},nf);

            % Symmetry enforcement
            if j ~= i
                model_temp{j, i} = model_temp{i, j};
                nu_temp{j, i} = nu_temp{i, j};
                c_temp{j, i} = c_temp{i, j};
            end
        end

        % Update directional asymmetry vectors if specified
        % (works only to linear variation)
        if isfield(param, 'dir')
            dir_values = computeDirectionalAsymmetry(param.dir{i}, tau(nf), cas, nbfam);
            if ~isempty(dir_values)
                vx(nf,i) = dir_values(1);
                if cas >= 2, vy(nf,i) = dir_values(2); end
                if cas >= 3, vz(nf,i) = dir_values(3); end
            end
        end
    end
    % Compute the spectral matrix and perform the eigenvalue–eigenvector decomposition
    % Negative eigenvalues are automatically replaced by 0, and the spectral
    % matrix is recomputed accordingly. Pseudo-admissibility needs to be
    % validated using TASC3D.m
    G{nf} = ComputeS(x0, x0r, model_temp, c_temp, nu_temp, Nx, Ny, Nz);

end
% Compute the spectral matrix and perform eigen-decomposition
% of the covariance model estimated from the data
G_mean = ComputeS(x0, x0r, model, c, nu, Nx, Ny, Nz);

% Compute the random number (randomly generated if no conditioning,
% carefully selected if conditioning is required). 'cexp' is the covariance
% matrix of the multivariate Gaussian random field generated to extract U.
[U, cexp] = GenerateU(G_mean, model, c, nu, param, nbsimul, nx, ny, nz, vx, vy, vz, tau, cas);

% Initialized Ubest
UBest = U;
% Clear some memory space
clear model_temp c_temp nu_temp dx dy dz x0 x0r Nx Ny Nz nf Nx2 Ny2 Nz2 G_mean model_local c_local nu_local

% Perform nbsimul simulations (the program works with the parfor loop)
for ns = 1:nbsimul

    % Initialize local variables for this iteration
    U_local = U{ns};                   % Local copy of U
    ErrBest_local = NaN(ittmax, 1);        % Local error tracking
    err = 10;                          % Initialize error
    
    % Compute the multivariate non-Gaussian random fields
    Z_local = ComputeField(G, U_local, cexp{ns}, param, tau, nbfam, nx, ny, nz, vx, vy, vz);

    ZBest_local = Z_local;
    UBest_local = U_local;

    % If required, perform Monte-Carlo optimization-based conditioning
    if ~isempty(param.HD)
        count = 1;       % Iteration counter 
        radius = 1;     % radius search for the perturbed point (you can adjust this)
        
        ErrBest_local(count) = ComputeError(Z_local, param.HD, nx, ny, cas);
        
        while err > 0.01 && count < ittmax
                [Z_opt, U_opt, err_opt] = OptimizeField(Z_local, G, U_local, cexp{ns}, tau, nbfam, param, ceil(radius), nx, ny, nz, vx, vy, vz, cas);
                count = count + 1;
                
                % Accept everything but keep the best solution
                if err_opt < err
                    UBest_local = U_opt;
                    ZBest_local = Z_opt;
                    err = err_opt;
                end                
                U_local = U_opt;
                Z_local = Z_opt;

                % Log error for this iteration
                ErrBest_local(count) = err;
        end
    end

    % Assign outputs to the sliced variables
    ZBest(:, ns, :) = ZBest_local;    % Assign ZBest for this simulation
    UBest{ns} = UBest_local;              % Assign the best U
    ErrBest(:, ns) = ErrBest_local;   % Assign error history
end

%Clear some memory space
clear U G i j cas nbsimul ns vx vy vz dir_values tau

% Post-conditioning by simple cokriging
if ~isempty(param.HD)
    ZBest = postcond(ZBest, param, model, c, nu, nx, ny, nz);
end

end


%% Compute the field
function Z = ComputeField(G, U, cexp, param, tau, nbfam, nx, ny, nz, vx, vy, vz)
% ComputeField: Simulates Gaussian fields and transforms them into non-Gaussian fields.
%
% Inputs:
%   G     - Cell array containing spectral density matrices (nbfam x nvar x nvar).
%   U     - Cell array of Gaussian white noise fields (nvar x size of field).
%   tau   - Threshold values for transformations (nbfam x 1).
%   nbfam - Number of field families.
%   nx, ny, nz - Grid dimensions.
%   vx, vy, vz - Shifts in the grid dimensions (nvar x 1).
%
% Output:
%   Z     - Transformed fields with non-Gaussian dependencies (nx*ny*nz x nvar).

% Number of variables
nvar = length(U);

% Avoid numerical instability (if 0 or 1, its becomes -inf or inf)

tau(1) = 0.00000001; tau(end) = 0.99999999;

% Initialize fields
datasim = zeros(nx * ny * nz, nbfam, nvar);  % Simulated data for each variable

% Perform Fourier transform of the Gaussian white noise
for j = 1:nvar
    U{j} = fftn(U{j});  % Fourier transform of Gaussian noise
end

% Loop over field families (nbfam)
for nf = 1:nbfam
    % Initialize the spectral density product for each family
    GU = cell(nvar, 1);
    for i = 1:nvar
        GU{i} = 0;  % Initialize to zero for each variable
    end

    % Compute the spectral density for each family and variable
    for j = 1:nvar
        for i = 1:nvar
            GU{i} = GU{i} + G{nf}{i, j} .* U{j};  % Spectral density multiplication
        end
    end

    % Perform the inverse Fourier transform to get random fields in real space
    for j = 1:nvar
        data = real(ifftn(GU{j}));  % Inverse Fourier transform
        datasim(:, nf, j) = reshape(data((ceil(nx/2):ceil(3*nx/2)-1)+vx(nf,j),(ceil(ny/2):ceil(3*ny/2)-1)+vy(nf,j),(ceil(nz/2):ceil(3*nz/2)-1)+vz(nf,j)), [], 1); % Reshape and store the simulated data
    end

end

% Transform to Non-Gaussian fields
if nbfam == 1
    Z = datasim;
else
    Z = LinearInterp(datasim, norminv(tau), param.Below);
end

% Ensure Gaussian distribution
Z = anamor_multi(Z, nvar);

% Reconstruct the corelation, i.e., normalization (avoid artefact from the truncation)
if nbfam>1 && nvar >1  
    Z = normalized_corr(Z, cexp);
elseif nbfam>1 && nvar == 1
    Z = Z / mean(std(Z));
end

% If Gaussian simulation
Z(:,param.Gaus)=datasim(:,ceil((nbfam+1)/2),param.Gaus);


end

%% Generate the multivariate white noise fields
function [U, cexp] = GenerateU(G, model, c, nu, param, nbsimul, nx, ny, nz, vx, vy, vz, tau, cas)
% GenerateU: Generates random fields U.
%
% This function generates random fields U for a simulation based on a given model,
% covariance matrix, and shape parameter. Depending on the conditioning setup,
% U is either generated from a conditional Gaussian realization (using post-conditioppning by cokriging)
% or from a random normal distribution. The function can handle different
% grid sizes and simulation configurations for multiple variables.
%
% Inputs:
%   G - Spectral decomposition (given as a matrix or list of matrices).
%   model - nvar x nvar cell array of structures defining the covariance models.
%             Each cell {i, j} defines the cross-covariance between variables i and j:
%             - 1D or isotropic: [model type, range].
%             - 2D anisotropic: [model type, range1, range2, rotx].
%             - 3D anisotropic: [model type, range1, range2, range3, rotx, roty, rotz].
%   c  - nvar x nvar cell array of sills defining the covariance models.
%   nu - nvar x nvar cell array of shape parameters for the covariance models.
%   nvar - Number of variables to be simulated.
%   nx, ny, nz - Dimensions of the simulation grid (x, y, z directions).
%   cas - Case parameter to define the dimensionality of the random field:
%         cas == 1 for 1D, cas == 2 for 2D, cas == 3 for 3D.
%
% Outputs:
%   U - A cell array of random fields (1D, 2D, or 3D depending on the case).
%

% Enlarge grid size to avoid aliasing
Nx = (nx + 1) * 2; Ny = (ny + 1) * 2; Nz = (nz + 1) * 2;

% Number of variables
nvar = size(model, 1);

% If conditioning is required (param.HD is provided), generate U using conditional Gaussian realizations.
if ~isempty(param.HD)

    U = cell(nbsimul,1); cexp = cell(nbsimul,1);
    % We generate unconditional simulation using GFFTMA.
    datasim = zeros(length(reshape(G{1},[],1)),nbsimul, nvar);
    for ns = 1 :nbsimul
        % Generate Fourier transforms of random Gaussian fields
        GU = cell(nvar, 1);
        for j = 1:nvar
            GU{j} = 0;
        end
        for j = 1:nvar
            if cas==1
                U{ns}{j,1} = randn(Nx,1);  % Gausssian white noise
            elseif cas ==2
                U{ns}{j,1} = randn(Nx, Ny);  % Gausssian white noise
            elseif cas ==3
                U{ns}{j,1} = randn(Nx, Ny, Nz);  % Gausssian white noise
            end
            for i = 1:nvar
                GU{i} = GU{i} + G{i, j} .* fftn(U{ns}{j});  % Combine covariance with Gausssian white noise (spectral domain)
            end
        end
        for j = 1:nvar
            datasimTemp = real(ifftn(GU{j}));  % Simulated data in spatial domain
            datasim(:,ns, j) = reshape(datasimTemp, [], 1);  % Flatten to column vector
        end
        % Direct- and Cross - Covariance of the experimental mean familly
        cexp{ns} = cov(squeeze(datasim(:, ns, :)));
    end
    clear datasimTemp GU
    %Post-conditioning using simple cokriging
    %Reposition the conditioning data at the middle
    if cas == 1
        [loc, ~] = computeLocationsAndIndices(param, tau, nx, 0, 0, vx, [], [], Nx, [], cas);
        for i = 1:nvar
            param.HD{i}(:, 1) = loc{i}(:, 1); % Update param.HD for x-coordinate
        end
        datasim = postcond(datasim, param, model, c, nu, Nx, 1, 1);

    elseif cas == 2
        [loc, ~] = computeLocationsAndIndices(param, tau, nx, ny, 0, vx, vy, [], Nx, Ny, cas);
        for i = 1:nvar
            param.HD{i}(:, 1) = loc{i}(:, 1); % Update param.HD for x-coordinate
            param.HD{i}(:, 2) = loc{i}(:, 2); % Update param.HD for y-coordinate
        end
        datasim = postcond(datasim, param, model, c, nu, Nx, Ny, 1);

    elseif cas == 3
        [loc, ~] = computeLocationsAndIndices(param, tau, nx, ny, nz, vx, vy, vz, Nx, Ny, cas);
        for i = 1:nvar
            param.HD{i}(:, 1) = loc{i}(:, 1); % Update param.HD for x-coordinate
            param.HD{i}(:, 2) = loc{i}(:, 2); % Update param.HD for y-coordinate
            param.HD{i}(:, 3) = loc{i}(:, 3); % Update param.HD for z-coordinate
        end
        datasim = postcond(datasim, param, model, c, nu, Nx, Ny, Nz);
    end
    
    GU =cell(nvar,1); U =cell(nbsimul,1);
    % Reconstruct the Gaussian random noise
    for ns = 1 : nbsimul
        if cas == 1
            for j = 1 : nvar
                GU{j} = fftn(datasim(:,ns,j));
            end
        elseif cas == 2
            for j = 1 : nvar
                GU{j} = fftn(reshape(datasim(:,ns,j),[Nx, Ny]));
            end
        elseif cas == 3
            for j = 1 : nvar
                GU{j} = fftn(reshape(datasim(:,ns,j),[Nx, Ny, Nz]));
            end
        end
        U{ns} = reconstruct_U(G, GU);
    end

    % If conditioning is not required, generate U from a normal distribution.
elseif isempty(param.HD)

    U = cell(nbsimul,1); cexp = cell(nbsimul,1);
    % We generate unconditional simulation using GFFTMA.
    datasim = zeros(length(reshape(G{1},[],1)),nbsimul, nvar);

    for ns = 1 :nbsimul
        % Generate Fourier transforms of random Gaussian fields
        GU = cell(nvar, 1);
        for j = 1:nvar
            GU{j} = 0;
        end
        for j = 1:nvar
            if cas==1
                U{ns}{j,1} = randn(Nx,1);  % Gausssian white noise
            elseif cas ==2
                U{ns}{j,1} = randn(Nx, Ny);  % Gausssian white noise
            elseif cas ==3
                U{ns}{j,1} = randn(Nx, Ny, Nz);  % Gausssian white noise
            end
            for i = 1:nvar
                GU{i} = GU{i} + G{i, j} .* fftn(U{ns}{j});  % Combine covariance with Gausssian white noise (spectral domain)
            end
        end
        for j = 1:nvar
            datasimTemp = real(ifftn(GU{j}));  % Simulated data in spatial domain
            datasim(:,ns, j) = reshape(datasimTemp, [], 1);  % Flatten to column vector
        end
        % Direct- and Cross - Covariance of the experimental mean familly
        cexp{ns} = cov(squeeze(datasim(:, ns, :)));
    end
else
    warning('param.HD is not defined correctly.');
end

end

%% Post-conditioning by simple cokriging
function datasim = postcond(datasim, param, model, c, nu, nx, ny, nz)
% postcond: Performs post-conditioning by simple cokriging.
%
% Inputs:
%   datasim - Simulated data matrix (N x nVar).
%   param   - A structure containing conditioning parameters, including param.HD.
%   model   - Covariance model definitions.
%   c       - Covariance model sills.
%   nu      - Covariance model shape parameters.
%   nx, ny, nz - Grid dimensions.
%
% Output:
%   paramStruct - A structure containing all necessary information for post-conditioning.

% Determine simulation dimensionality and grid points
if nz == 1 && ny == 1 && nx > 1
    x0 = (1:nx)';  % 1D grid
    d = 1;
elseif nz == 1 && ny > 1 && nx > 1
    x0 = grille2(1, nx, 1, 1, ny, 1);  % 2D grid
    d = 2;
elseif nz > 1 && ny > 1 && nx > 1
    x0 = grille3(1, nx, 1, 1, ny, 1, 1, nz, 1);  % 3D grid
    d = 3;
else
    error('Check input parameters: Incorrect grid dimensions.');
end

nvar = size(model, 1);  % Number of variables
x0grid = cell(nvar, 1); x = cell(nvar, 1); xgrid = cell(nvar, 1);
idx = cell(nvar, 1); xsim = cell(nvar, 1); x0sim = cell(nvar, 1);

% Loop over the number of variables
for i = 1:nvar
    x{i, 1} = param.HD{i}(~isnan(param.HD{i}(:, d + 1)), d + 1);
    xgrid{i, 1} = param.HD{i}(~isnan(param.HD{i}(:, d + 1)), 1:2);
    idx{i} = (xgrid{i, 1}(:, 2) - 1) * nx + xgrid{i, 1}(:, 1);
    xsim{i, 1} = datasim(idx{i},:, i);
    x0sim{i, 1} = datasim(:,:, i);
end

% Precompute covariance matrices of observed points
K = ComputeS_Cokri(xgrid, xgrid, model, c, nu, 1);

% Extract precomputed values
nbsimul = size(x0sim{1},2);

% Vectorize inputs
x = cell2mat(x);
xsim = cell2mat(xsim);

% Compute residuals
xres = x - xsim;

% Solve the kriging system (dual form)
bres = K \ xres;

% Iteratively cokriged the residuals by batch
batchSize = 10000; % Define the size of each batch
M = size(x0,1); % Total number of simulated point
numBatches = ceil(M/ batchSize); % Calculate the total number of batches

%Reinitialized datasim
datasim = zeros(M,nbsimul,nvar);
for i = 1 : nvar
    for batchIdx = 1:numBatches
        % Calculate the start and end column indices for the current batch
        startCol = (batchIdx - 1) * batchSize + 1;
        endCol = min(batchIdx * batchSize, M);        
        for j=1:nvar
            if i == j
                x0grid{j}=x0(startCol:endCol,:);
            else
                x0grid{j}=[];
            end
        end
        % Covariance matrices of observed points with grid point
        k = ComputeS_Cokri(xgrid, x0grid, model, c, nu, 0);

        % Compute cokriged residuals
        xres_krig = (bres' * k)';

        % Combine unconditional simulations and cokriged residuals
        datasim(startCol:endCol,:,i) = xres_krig + x0sim{i}(startCol:endCol,:);
    end
end

end

%% Helper to reconstruct the Gaussian white noise from a conditioned random field and a specified covariance matrix.
function U = reconstruct_U(G, GU)
% G is a nvar x nvar cell array of covariance matrices (M entries each)
% GU is a nvar x 1 cell array of cross-covariance products (M entries each)
% cas indicates dimensionality (1 for 1D, 2 for 2D, 3 for 3D)

% Number of variables
nvar = size(G, 1);

% Determine the total number of points (M)
M = numel(G{1, 1});  % Total number of spatial points

% Convert cell array G into a 3D array: nvar x nvar x M
G_array = zeros(nvar, nvar, M);
for i = 1:nvar
    for j = 1:nvar
        G_array(i, j, :) = reshape(G{i, j}, 1, 1, M);  % Reshape to match dimensions
    end
end

% Convert GU into a 2D array: nvar x M
GU_array = zeros(nvar, M);
for i = 1:nvar
    GU_array(i, :) = reshape(GU{i}, 1, M);  % Reshape to match dimensions
end

% Batch inversion of G using pageinv
G_inv_array = pageinv(G_array);  % G_inv_array is nvar x nvar x M
G_inv_array(G_inv_array==inf)=0; % Correct for inversion of matrix of zeros

% Batch matrix multiplication of G_inv_array and GU_array
U_array = pagemtimes(G_inv_array, reshape(GU_array, [nvar, 1, M]));  % U_array is nvar x 1 x M

% Reshape U_array back to cell format
U = cell(nvar, 1);
for i = 1:nvar
    U{i} = reshape(U_array(i, 1, :), size(G{1, 1}));  % Reshape to original spatial dimensions
end

% Perform inverse FFT for all variables
for i = 1:nvar
    U{i} = ifftn(U{i});
end
end

%% Unified function to update parameters
function updatedValue = updateParameters(currentValue, param, tau_nf, i, j, nbfam, fieldMapping,nf)
updatedValue = currentValue;
for k = 1:size(fieldMapping, 1)
    field = fieldMapping{k, 1};
    idx = fieldMapping{k, 2};
    if isfield(param, field)
        fieldParam = param.(field){i, j};
        if length(fieldParam) == 2
            newValue = (fieldParam(2) - fieldParam(1)) * tau_nf + fieldParam(1);
        elseif length(fieldParam) == nbfam
            newValue = fieldParam(nf);
        else
            continue;
        end

        if ~isempty(idx)
            updatedValue(:,idx) = newValue; % Update specific index for models
        else
            updatedValue = newValue; % Update entire value for scalar parameters
        end
    end
end
end

%% Helper function to compute directional asymmetry
function dir_values = computeDirectionalAsymmetry(dir_param, tau_nf, cas, nbfam)
if size(dir_param, 2) == 2
    dir_values = round((dir_param(:, 2) - dir_param(:, 1)) * tau_nf + dir_param(:, 1));
elseif size(dir_param, 2) == nbfam
    dir_values = round(dir_param(:, nf));
elseif size(dir_param, 2) == 1
    dir_values = round(dir_param(:));
else
    dir_values = [];
end
dir_values = dir_values(1:min(cas, length(dir_values))); % Adjust for the number of cases
end

%% Postprocessing function to ensure the reproduction of the covariance matrix.
function data = normalized_corr(data,c)
% This function computes the average 2x2 correlation matrix across all time steps
% Input:
%   data - A 3D matrix of size N x M x 2, where N is the number of data points,
%          M is the number of time steps, and 2 corresponds to two variables.
% Output:
%   avg_corr_matrix - The average 2x2 correlation matrix over all time steps.

% Extract data for the current time step
data_at_time = squeeze(data);

% Compute the correlation matrix for this time step
corr_matrix = cov(data_at_time);

% Choleski decomposition
L_star = chol(corr_matrix,'lower');
L =chol(c,'lower');
LL = L /L_star;

% Perform the matrix multiplication for each 2D slice
data = ( LL * squeeze(data)' )';

end

%% %% Function to ensure that the cokriging matrix is semi-definite positive. Negative eigenvalues are replaced with 0.
function C = ComputeS_Cokri(x0, x0r, model, c, nu, idx)
% ComputeS: Function to compute the spectral densities of multivariate covariance matrices
% and perform Cholesky decomposition for spectral simulation.
%
% Syntax:
%   G = ComputeS_Cokri(x0, x0r, model)
%
% Description:
%   This function calculates the spectral decomposition of multivariate covariance matrices
%   for a given set of input parameters. It uses the covariance model to compute
%   Fourier-transformed covariance matrices and applies Cholesky decomposition
%   at all frequencies for spectral simulation.
%
% Inputs:
%   x0     : [n x d] matrix of input coordinates.
%   x0r    : [1 x d] matrix of the coordinates of the center grid.
%   model  : nvar x nvar cell array of structures defining the covariance models.
%             - Each cell {i, j} contains structures for cross-covariance
%               between variables i and j.
%             - Isotropic case: [model type, range, shape, sill].
%             - Anisotropic case: [model type, range1, range2, range3, rotx, roty, rotz, shape, sill].
%   nn     : Integer, number of realizations.
%
% Outputs:
%   G      : [nvar x nvar] cell array containing the Cholesky-decomposed
%            spectral matrices at all frequencies.
%
% Notes:
%   - Supports 1D, 2D, or 3D cases based on grid dimensions.
%   - Uses FFT to compute the Fourier transform of the covariance matrices.
%   - Handles isotropic and anisotropic covariance models.
%
% Author:
%   D. Marcotte, May 2015.
%   D. Lauzon , November 2024. (Modified from Marcotte, 2015)
%

% Number of variables
nvar = size(model, 1);

% Initialize covariance matrix
C = cell(nvar, nvar);

% Perform Cholesky decomposition
if idx == 1
    for i = 1:nvar
        for j = 1:i
            C{i, j} = covar(x0{i}, x0r{j}, model{i, j}, c{i, j}, nu{i, j});
            C{j, i} = C{i, j}'; % Symmetry
        end
    end
    C = cell2mat(C);
else
    for i = 1:nvar
        for j = 1:nvar
            if isempty(x0r{j})
                C{i, j}= [];
            else
                C{i, j} = covar(x0{i}, x0r{j}, model{i, j}, c{i, j}, nu{i, j});
            end
        end
    end
    C = cell2mat(C);
end
[v , d] = eig(C);
C = v * diag(max(diag(d), 0)) * v';
end

function G = ComputeS(x0, x0r, model, c, nu, Nx, Ny, Nz)
% ComputeS: Function to compute the spectral densities of multivariate covariance matrices
% and perform Cholesky decomposition for spectral simulation.
%
% Syntax:
%   G = ComputeS(x0, x0r, model, Nx, Ny, Nz)
%
% Description:
%   This function calculates the spectral decomposition of multivariate covariance matrices
%   for a given set of input parameters. It uses the covariance model to compute
%   Fourier-transformed covariance matrices and applies Cholesky decomposition
%   at all frequencies for spectral simulation.
%
% Inputs:
%   x0     : [n x d] matrix of input coordinates.
%   x0r    : [1 x d] matrix of the coordinates of the center grid.
%   model  : nvar x nvar cell array of structures defining the covariance models.
%             - Each cell {i, j} contains structures for cross-covariance
%               between variables i and j.
%             - Isotropic case: [model type, range, shape, sill].
%             - Anisotropic case: [model type, range1, range2, range3, rotx, roty, rotz, shape, sill].
%   nn     : Integer, number of realizations.
%   Nx     : Integer, number of grid points along the x-axis (for 1D cases).
%   Ny     : Integer, number of grid points along the y-axis (for 2D cases).
%   Nz     : Integer, number of grid points along the z-axis (for 3D cases).
%
% Outputs:
%   G      : [nvar x nvar] cell array containing the Cholesky-decomposed
%            spectral matrices at all frequencies.
%
% Notes:
%   - Supports 1D, 2D, or 3D cases based on grid dimensions.
%   - Uses FFT to compute the Fourier transform of the covariance matrices.
%   - Handles isotropic and anisotropic covariance models.
%
% Author:
%   D. Marcotte, May 2015.
%   D. Lauzon , November 2024. (Modified from Marcotte, 2015)
%

% Compute the case
cas = size(x0,2);

% Number of variables
nvar = size(model, 1);

% Initialize covariance matrix
S = cell(nvar, nvar);

for i = 1:nvar
    for j = 1:i
        cc = covar(x0, x0r, model{i, j}, c{i, j}, nu{i, j});

        % Handle different cases
        switch cas
            case 1
                cc = fftshift(cc);
            case 2
                cc = fftshift(reshape(cc, Nx + 1, Ny + 1));
            case 3
                cc = fftshift(reshape(cc, Nx + 1, Ny + 1, Nz + 1));
        end

        % Fourier transform of covariances
        S{i, j} = real(fftn(cc));
        S{j, i} = S{i, j}; % Symmetry
    end
end

% Perform Cholesky decomposition
G = chol_dec(S);

end

function G = chol_dec(S)
% chol_dec: Helper function for Cholesky decomposition of spectral matrices.
%
% Syntax:
%   G = chol_dec(S)
%
% Description:
%   Performs Cholesky decomposition of spectral matrices for use in
%   spectral simulation. This function handles matrices stored in a
%   cell array format and processes only the significant frequencies.
%
% Inputs:
%   S      : [nvar x nvar] cell array of spectral covariance matrices.
%            - Each cell contains a matrix of spectral coefficients
%              corresponding to a given variable pair.
%
% Outputs:
%   G      : [nvar x nvar] cell array of decomposed spectral matrices.
%            - Each cell contains a matrix of Cholesky decomposition
%              results at all frequencies.
%
% Notes:
%   - Applies eigenvalue decomposition to ensure non-negative eigenvalues.
%   - Processes only significant frequencies to improve computational efficiency.
%
% Author:
%   D. Marcotte, May 2015.
%   D. Lauzon , November 2024. Vectorized the code.
%
% Perform Cholesky decomposition simultaneously at all frequencies
% S is a cell array (nvar x nvar), each cell containing a matrix of spectrum coefficients

nvar = size(S, 1);
[nx, ny] = size(S{1, 1});
total_points = nx * ny;

% Identify significant indices
id = zeros(nx, ny);
for i = 1:nvar
    id = id | (S{i, i} > max(S{i, i}(:)) * 1e-6); % Logical OR for significant values
end
id = find(id); % Linear indices of significant points
num_points = length(id);

% Preallocate 3D array for significant indices: nvar x nvar x num_points
m = zeros(nvar, nvar, num_points);

% Populate m matrices for significant indices
for i1 = 1:nvar
    for i2 = i1:nvar
        values = reshape(S{i1, i2}, total_points, 1);
        m(i1, i2, :) = values(id);
        if i1 ~= i2
            m(i2, i1, :) = values(id);
        end
    end
end

% Perform eigenvalue decomposition in batch
[V, D] = pageeig(m); % Batch eigenvalue decomposition (nvar x nvar x num_points)

% Construct G using eigen decomposition results
D = max(D, 0); % Ensure non-negative eigenvalues
G_matrices = pagemtimes(V, pagemtimes(sqrt(D), pagetranspose(V))); % G = V * sqrt(D) * V'

% Reassign results back to G
G = cell(nvar, nvar);
for i1 = 1:nvar
    for i2 = 1:nvar
        temp = zeros(total_points, 1);
        temp(id) = squeeze(G_matrices(i1, i2, :));
        G{i1, i2} = reshape(temp, nx, ny);
    end
end
end

%% Monte Carlo optimization for conditioning (not perfect).
function [Z_opt, U_opt, err_opt] = OptimizeField( Z, G, U, cexp, tau, nbfam, param, radius, nx, ny, nz, vx, vy, vz, cas)
    % OptimizeField: Optimizes the field using perturbation and error minimization to find the optimal perturbation t.
    % This function uses fminsearch to minimize the error and return the optimized field.
    %
    % Inputs:
    %   G       - Cell array of spectral densities for the random fields.
    %   U       - Cell array of Gaussian white noise fields.
    %   tau     - Threshold levels for the transformation.
    %   nbfam   - Number of random field families.
    %   paramHD - Cell array containing target parameters (indices and target values).
    %   nx, ny, nz - Dimensions of the field.
    %   vx, vy, vz - Shifts to apply during field computation.
    %
    % Outputs:
    %   Z_opt   - Optimized field based on the perturbation t.
    %   U_opt   - Optimized Gaussian white noise fields.
    %   err_opt - Error (RMSE) associated with the optimized field.
    %   t_opt   - Optimal perturbation angle in radians that minimizes the error.
    
    % Number of variables (fields)
    nvar = size(U, 1);
    
    % Generated the fusion fields
    u = cell(size(U));
    % Extract indices and target values
    for j = 1:nvar
        if cas ==1
            indices = param.HD{j}(:, 1);
        elseif cas==2
            indices = (param.HD{j}(:, 2)-1)*nx + param.HD{j}(:, 1);
        elseif cas ==3
            indices = (param.HD{j}(:, 3)-1)*nx*ny + (param.HD{j}(:, 2)-1)*nx + param{j}(:, 1);
        else
            warning('Error in dimension or dimension not support')
        end

        u{j} = param.HD{j}(:,end) - Z(indices,j);
    end
    [Z_opt, U_opt, err_opt] = PerturbeField( G, U, u, cexp, tau, nbfam, param, radius, nx, ny, nz, vx, vy, vz, cas);

end

%% Applied the perturbation of the Monte Carlo optimization for conditioning (not perfect).
function [Z, U, err] = PerturbeField( G, U, u, cexp, tau, nbfam, param, radius, nx, ny, nz, vx, vy, vz, cas)
    % PerturbeFieldWithT: Perturbs the Gaussian white noise with a given angle t
    % and computes the resulting random field and its associated error.
    %
    % Inputs:
    %   G       - Cell array of spectral densities for the random fields.
    %   U       - Cell array of Gaussian white noise fields.
    %   tau     - Threshold levels for the transformation.
    %   nbfam   - Number of random field families.
    %   paramHD - Cell array of target parameters (indices and values).
    %   nx, ny, nz - Dimensions of the field.
    %   vx, vy, vz - Shifts to apply during field computation.
    %   t       - Perturbation angle in radians.
    %
    % Outputs:
    %   Z       - Simulated random field.
    %   U       - Perturbed Gaussian white noise fields.
    %   err     - Root Mean Squared Error (RMSE) between simulated and target values.
    
    % Number of variables (fields)
    nvar = size(U, 1);

    % Enlarge grid size to avoir aliasing
    Nx=(nx+1)*2; Ny=(ny+1)*2;
    
    % Validate input dimensions
    if ~iscell(U) || ~iscell(G)
        error('Inputs U and G must be cell arrays.');
    end
    
    if length(G) < nbfam
        error('Number of families (nbfam) exceeds the size of G.');
    end
    
    [loc, ~] = computeLocationsAndIndices(param, tau, nx, ny, nz, vx, vy, vz, Nx, Ny, cas);
    
    for i = 1:nvar
        pdf = reshape(radial_decay(size(U{i}), loc{i}, u{i}.*abs(rand(size(u{i})))*0.75, radius),[],1);
        U{i}= U{i} + reshape(pdf,[Nx Ny]);
    end

    % Compute the random field Z
    Z = ComputeField(G, U, cexp, param, tau, nbfam, nx, ny, nz, vx, vy, vz);
    
    % Compute the error (RMSE) compared to target values
    err = ComputeError(Z, param.HD, nx, ny, cas);
end

%% Function to compute the misfit between real values and simulated value
function err = ComputeError(Z, param, nx, ny, cas)
    % ComputeError: Computes the Root Mean Squared Error (RMSE) based on provided indices and target values.
    %
    % Inputs:
    %   Z       - Simulated field (matrix of size [n_points, n_vars]).
    %   paramHD - Cell array where each cell contains a matrix with two columns:
    %             - First column: Indices for the points in Z.
    %             - Second column: Target values to compare with.
    %
    % Output:
    %   err     - Computed RMSE as a scalar.
    
    % Validate inputs
    if ~iscell(param)
        error('paramHD must be a cell array.');
    end
    
    % Number of variables
    nvar = length(param);
    
    % Initialize the error
    err = 0;
    
    % Compute the RMSE for each variable
    for j = 1:nvar
        % Extract indices and target values
        if cas ==1
            indices = param{j}(:, 1);
        elseif cas==2
            indices = (param{j}(:, 2)-1)*nx + param{j}(:, 1);
        elseif cas ==3
            indices = (param{j}(:, 3)-1)*nx*ny + (param{j}(:, 2)-1)*nx + param{j}(:, 1);
        else
            warning('Error in dimension or dimension not support')
        end
    
        % Validate indices
        if any(indices < 1) || any(indices > size(Z, 1))
            error('Indices in paramHD{%d} are out of bounds for the matrix Z.', j);
        end
        
        % Compute the RMSE for the current variable
        err = err + sqrt(mean((Z(indices, j) - param{j}(:, end)).^2));
    end
    
    % Normalize by the number of variables
    err = err / nvar;
end

%% Function to transform the multivariate Gaussian random field set to a non-Gaussian field
function x_crossing = LinearInterp(datasim, tau, param)
% LinearInterp - Calculates the crossing points between two linear segments
% based on the input matrix `datasim` and threshold vector `tau`.
%
% Syntax:
%   x_crossing = LinearInterp(datasim, tau)
%
% Inputs:
%   - datasim: n x m matrix, where each row contains a smooth function from 
%     element 1 to element m.
%   - tau: 1 x m vector representing the threshold values that the elements
%     of datasim should be compared to.
%
% Outputs:
%   sing: n x 1 vector containing the interpolated crossing points 
%     where the curve crosses the threshold.

% Example Data- x_cros
[n, m, nvar] = size(datasim);

% Initialize output with NaN
x_crossing = NaN(n, nvar);

for j = 1 : nvar
    datasim_now = datasim(:, :, j);

    if param(j) == 1  % Bas → haut (montée)
        indicator_matrix = (datasim_now(:, 1:end-1) >= tau(1:end-1)) & ...
            (datasim_now(:, 2:end)   <  tau(2:end));

        [rows, cols] = find(indicator_matrix);
        first_crossing = accumarray(rows, cols, [n, 1], @max, 1);
        index_i_plus_1 = min(max(first_crossing+1, 1), m);

    elseif param(j) == 2  % Haut → bas (descente)
        datasim_now = datasim_now(:,end:-1:1);
        indicator_matrix = (datasim_now(:, 1:end-1) >= tau(1:end-1)) & ...
            (datasim_now(:, 2:end)   <  tau(2:end));

        [rows, cols] = find(indicator_matrix);
        first_crossing = accumarray(rows, cols, [n, 1], @min, 1);
        index_i_plus_1 = min(max(first_crossing+1, 1), m);
    end

    valid_crossing = first_crossing > 0;
    valid_rows = find(valid_crossing);

    A1 = datasim_now(sub2ind(size(datasim_now), valid_rows, first_crossing(valid_crossing)));
    A2 = datasim_now(sub2ind(size(datasim_now), valid_rows, index_i_plus_1(valid_crossing)));
    tau1 = tau(first_crossing(valid_crossing))';
    tau2 = tau(index_i_plus_1(valid_crossing))';

    a = (A1 - A2) ./ (tau1 - tau2);
    b = A1 - a .* tau1;
    x_crossing(valid_crossing, j) = b ./ (1 - a);
end

end

%% Radial window for Monte Carlo optimization
function map = radial_decay(grid_size, locations, values, radius)
    % Generalized radial decay map for 1D, 2D, and 3D grids.
    % Inputs:
    %   grid_size: Array specifying the size of the grid ([nx], [nx, ny], or [nx, ny, nz])
    %   locations: Matrix of size [n_locations, n_dims] (coordinates of the points)
    %   values: Vector of size [n_locations] (values at each location)
    %   radius: Radius of influence for the decay
    % Output:
    %   map: Grid with decayed values

    % Initialize the map
    map = zeros(grid_size);

    % Dimensionality of the problem
    n_dims = length(grid_size);

    % Avoid potential numerical instability
    values(values==0)=10^-8;

    % Loop through each location
    for k = 1:size(locations, 1)
        % Get the location and value
        loc = locations(k, :);
        value = values(k);

        % Define the bounding box around the location within the radius
        bounds_min = max(1, loc - radius);
        bounds_max = min(grid_size, loc + radius);

        % Generate relative coordinates within the patch
        switch n_dims
            case 1
                % 1D case
                x = bounds_min(1):bounds_max(1);
                dx = x - loc(1);
                dist = abs(dx);
                mask = dist < radius;
                decay = value * exp(-dist(mask).^2 / (0.1 * radius^2)) ;
                decay = decay/sum(decay)*value;
                map(x(mask)) = map(x(mask)) + decay;

            case 2
                % 2D case
                [dy, dx] = ndgrid((bounds_min(1):bounds_max(1)) - loc(1), ...
                                  (bounds_min(2):bounds_max(2)) - loc(2));
                dist = sqrt(dx.^2 + dy.^2);
                mask = dist < radius;
                decay = value * exp(-dist(mask).^2 / (0.1 * radius^2));
                decay = decay/sum(decay)*value;
                temp = map(bounds_min(1):bounds_max(1), bounds_min(2):bounds_max(2));
                temp(mask) = temp(mask) + decay;
                map(bounds_min(1):bounds_max(1), bounds_min(2):bounds_max(2)) = temp;

            case 3
                % 3D case
                [dz, dy, dx] = ndgrid((bounds_min(1):bounds_max(1)) - loc(1), ...
                                      (bounds_min(2):bounds_max(2)) - loc(2), ...
                                      (bounds_min(3):bounds_max(3)) - loc(3));
                dist = sqrt(dx.^2 + dy.^2 + dz.^2);
                mask = dist < radius;
                decay = value * exp(-dist(mask).^2 / (0.1 * radius^2)) ;
                decay = decay/sum(decay)*value;
                temp = map(bounds_min(1):bounds_max(1), bounds_min(2):bounds_max(2), bounds_min(3):bounds_max(3));
                temp(mask) = temp(mask) + decay;
                map(bounds_min(1):bounds_max(1), bounds_min(2):bounds_max(2), bounds_min(3):bounds_max(3)) = temp;

            otherwise
                error('Unsupported dimensionality: %d. Only 1D, 2D, and 3D are supported.', n_dims);
        end
    end
end

%% Compute the index to map from enlarged fields to original fields. It also takes into account the advective vector.

function [loc, idx] = computeLocationsAndIndices(param, tau, nx, ny, nz, vx, vy, vz, Nx, Ny, cas)
% computeLocationsAndIndices - Computes the locations and indices for Gaussian noise perturbation.
%
% Syntax:
%   [loc, idx] = computeLocationsAndIndices(param, tau, nx, ny, nz, vx, vy, vz, Nx, Ny, cas)
%
% Inputs:
%   - param: A structure containing HD data for each variable.
%   - tau: Threshold value for determining tagVec using the Gaussian inverse cumulative distribution.
%   - nx, ny, nz: Dimensions of the domain.
%   - vx, vy, vz: Functions or arrays to perturb the locations based on tagVec.
%   - Nx, Ny: Grid dimensions (number of elements along x and y).
%   - cas: Case selector (1, 2, or 3).
%
% Outputs:
%   - loc: A cell array containing location coordinates for each variable.
%   - idx: A cell array containing computed indices for each variable.

    nvar = length(param.HD); % Number of variables
    loc = cell(nvar, 1);
    idx = cell(nvar, 1);

    for j = 1:nvar
        % Compute tagVec based on the threshold tau
        if isscalar(tau)
            tau = 0;
        end
        tagVec = sum(param.HD{j}(:, end) > norminv(tau), 2);
        
        % Initialize loc for the current variable
        switch cas
            case 1
                % Case 1: Only x-coordinate
                loc{j}(:, 1) = param.HD{j}(:, 1) + ceil(nx / 2) - 1 + vx(tagVec, j);
                idx{j} = loc{j}(:, 1);
            case 2
                % Case 2: x- and y-coordinates
                loc{j}(:, 1) = param.HD{j}(:, 1) + ceil(nx / 2) - 1 + vx(tagVec, j);
                loc{j}(:, 2) = param.HD{j}(:, 2) + ceil(ny / 2) - 1 + vy(tagVec, j);
                idx{j} = (loc{j}(:, 2) - 1) * Nx + loc{j}(:, 1);
            case 3
                % Case 3: x-, y-, and z-coordinates
                loc{j}(:, 1) = param.HD{j}(:, 1) + ceil(nx / 2) - 1 + vx(tagVec, j);
                loc{j}(:, 2) = param.HD{j}(:, 2) + ceil(ny / 2) - 1 + vy(tagVec, j);
                loc{j}(:, 3) = param.HD{j}(:, 3) + ceil(nz / 2) - 1 + vz(tagVec, j);
                idx{j} = (loc{j}(:, 3) - 1) * Nx * Ny + (loc{j}(:, 2) - 1) * Nx + loc{j}(:, 1);
            otherwise
                error('Invalid case selected. Choose cas = 1, 2, or 3.');
        end
    end
end