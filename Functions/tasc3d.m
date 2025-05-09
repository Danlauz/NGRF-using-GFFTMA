function [icode, ics] = tasc3d(model, c, nu, s)
% TASC3D Tests the admissibility of a non-LMC model.
%
% Syntax:
%   [icode, ics] = tasc3d(model, c, nu, s)
%
% Inputs:
%   model - nvar x nvar cell array, where model{i,j} = model{j,i}.
%           Each cell contains n_ij structures defined by:
%           - Isotropic case: [type, range]
%           - Anisotropic case: [type, range1, range2, range3, rotx, roty, rotz]
%               - Rotations are counterclockwise: first around z, then y, then x.
%           Parameters may be set to NaN for estimating admissible values.
%           Note: If more than two parameters are NaN, execution stops with a message.
%
%   c     - nvar x nvar cell array of sill values, where c{i,j} = c{j,i}.
%           Each cell contains n_ij structures defined by a single parameter.
%
%   nu    - nvar x nvar cell array of shape parameters, where nu{i,j} = nu{j,i}.
%           Each cell contains n_ij structures defined by a single parameter:
%           - Alpha for Cauchy models.
%           - Nu for K-Bessel models.
%
%   s     - Vector of frequencies for admissibility testing.
%           - For anisotropic cases, this defines a 3D frequency grid.
%           - Example: s = [0.005:0.005:0.01, 0.02:0.02:0.1, 0.2:0.2:7]';
%           - Note: Large vectors result in larger grids (e.g., 100 frequencies = 1M grid nodes).
%
% Outputs:
%   icode - Admissibility flag:
%           - 1 if the model is admissible.
%           - 0 if the model is not admissible.
%
%   ics   - Cauchy-Schwarz condition flag:
%           - 1 if the model satisfies the Cauchy-Schwarz condition.
%           - 0 otherwise.
%           - Always: icode <= ics.
%
% Notes:
%   - The program identifies two domains:
%       1. The Cauchy-Schwarz admissible domain (necessary condition).
%       2. The real admissible domain.
%
% Authors:
%   D. Marcotte, April 2015
%   D. Lauzon, November 2024


% Initialization of parameters
nvar = size(model, 1);  % Number of variables
an = -1;               % Default flag for anisotropic models
cn = -1;               % Default flag for isotropic models

% Determine if isotropic or anisotropic model
ncol = size(model{1,1}, 2);  
iso = (ncol < 7);      % Flag for isotropic case (true if fewer than 7 parameters)

% Set default frequencies for admissibility testing if not provided
if isempty(s)
    s = 10.^(-3 : 0.1 : 2)';
end

% Index for range iso : [itype range], aniso : [itype range1 range2 range3]
if iso
    ir = 2;   % Index for range in isotropic case
else
    ir = 4;   % Index for range in anisotropic case
end

% Find the largest range and the largest nugget effect across all models
for i = 1:nvar
    for j = i:nvar
        % Update 'an' with the maximum range value found
        an = max(an, max(max(model{i, j}(:, 2:ir))));         
        % Update 'cn' with the maximum nugget effect found
        cn = max(cn, max(c{i, j})); 
    end
end


% Normalize model and covariance parameters
for i = 1:nvar
    for j = i:nvar
        % Normalize the model parameters (excluding the first column)
        model{i, j}(:, 2:ir) = model{i, j}(:, 2:ir) / an;        
        % Ensure symmetry by copying the normalized values to the transposed position
        model{j, i} = model{i, j};        
        % Normalize the covariance values
        c{i, j} = c{i, j} / cn;        
        % Ensure symmetry by copying the normalized covariance to the transposed position
        c{j, i} = c{i, j};
    end
end

% Separate treatment for the nugget effect
% First, check if the nugget effect is admissible. If not, icode=0 and stop.
icode = check_nugget(model, c);
ics = [];
if icode == 0  % Stop if nugget effect is not admissible
    ics = nan;  % Set ics to NaN to indicate failure
    return;     % Exit the function early
end

% Eliminate the nugget effect from the spectral density test
for i = 1:nvar
    for j = 1:nvar
        % Identify components that are different from the nugget effect (first column > 1)
        id = model{i, j}(:, 1) > 1;        
        % Keep only the components that are not related to nugget effect
        model{i, j} = model{i, j}(id, :); c{i, j} = c{i, j}(id, :); nu{i, j} = nu{i, j}(id, :);
        % Ensure symmetry by copying the updated model values to the transposed position
        model{j, i} = model{i, j}; c{j, i} = c{i, j}; nu{j, i} = nu{i, j};
    end
end


% Test if the model satisfies the Cauchy-Schwarz condition
ics = test_cs(model, c, nu, iso);  % Test if the model satisfies the Cauchy-Schwarz condition
% If the Cauchy-Schwarz condition is satisfied
% Check if the model is admissible
% If the Cauchy-Schwarz condition is not satisfied, set icode to 0
if ics == 1
    icode = check_admiss(model,c, nu, iso, s);
else
    icode = 0;
end

function ics = test_cs(model, c, nu, iso)
% This function checks the Cauchy-Schwarz condition for a given model.
%
% Syntax:
%   ics = test_cs(model, c, nu, iso)
%
% Inputs:
%   - model: A nvar x nvar cell array representing the model for each pair (i,j).
%   - c: A nvar x nvar cell array representing the sill for each pair (i,j).
%   - nu: A nvar x nvar cell array representing the shape parameters for each pair (i,j).
%   - iso: A boolean flag indicating whether the model is isotropic (true) or anisotropic (false).
%
% Outputs:
%   - ics: 1 if the model satisfies the Cauchy-Schwarz condition, 0 otherwise.

% Define a vector of normalized distances where the CS condition is checked
h = [0, 0.001:0.001:0.01, 0.02:0.01:0.1, 0.2:0.1:1, 2:1:10, 15:5:50]'; 

% For isotropic models, the grid is one-dimensional; for anisotropic, it's 3D
if ~iso
    [h1, h2, h3] = ndgrid(h, h, h);  % Create 3D grid of distances for anisotropic case
    h = [h1(:), h2(:), h3(:)];        % Flatten the grid into a 3-column matrix
end

nvar = size(model, 1);  % Number of variables in the model
ics = 1;  % Assume the model meets the CS condition initially

% Loop through each pair of models (i, j)
for i = 1:nvar-1
    % Compute the term for model(i,i) - the auto-covariance term
    g11 = sum(c{i,i}) - covar(h, zeros(1, size(h, 2)), model{i,i}, c{i,i}, nu{i,i});
    
    for j = i+1:nvar
        % Compute the terms for model(j,j) and the cross-covariance term (i,j)
        g22 = sum(c{j,j}) - covar(h, zeros(1, size(h, 2)), model{j,j}, c{j,j}, nu{j,j});
        g12 = sum(c{i,j}) - covar(h, zeros(1, size(h, 2)), model{i,j}, c{i,j}, nu{i,j});
        
        % Compute the Cauchy-Schwarz discrepancy for the current pair (i,j)
        d = sqrt(g11 .* g22) - abs(g12);
        
        % If the condition is violated (i.e., d is negative), return 0
        if min(d) < -1e-09
            ics = 0;  % One of the models does not meet the CS condition
            return;
        end
    end
end

function icode = check_admiss(model, c, shape, iso, s)
% This function checks the admissibility of a multivariate covariance model
% by evaluating the spectral densities for the given models.
%
% Syntax:
%   icode = check_admiss(model, c, shape, iso, s)
%
% Inputs:
%   - model: A nvar x nvar cell array containing model parameters for each pair (i,j).
%   - c: A nvar x nvar cell array containing the sill (covariance) values for each pair.
%   - shape: A nvar x nvar cell array containing the shape parameter for each model.
%   - iso: A boolean flag indicating whether the model is isotropic (true) or anisotropic (false).
%   - s: A vector of frequency values at which to test admissibility.
%
% Outputs:
%   - icode: 1 if the model is admissible, 0 otherwise.

nvar = size(model, 1);
icode = 1;  % Assume the model is admissible initially

% Define spectral densities for various models (N/A implies not
% implemented, or no densities are admitted)
nugget = '0'; 
expon  = 'at/pi^2./(1+sa.^2).^2';
gaus   = 'at/(8*pi^1.5)*exp(-0.25*sa.^2)';
spher  = '0.75*at/pi./sa.^3.*besselj(1.5,sa/2).^2';
cubic  = '210*at./(pi^2*sa.^10).*(6*sa.*cos(sa/2)+(sa.^2-12).*sin(sa/2)).^2';
penta  = '27720*at./(pi^2*sa.^14).*((sa.^3-60*sa).*cos(sa/2)+(120-12*sa.^2).*sin(sa/2)).^2';
cauchy = 'at*sa.^(nu-1.5)/(pi^1.5*2^(nu+0.5)*gamma(nu)).*besselk(1.5-nu,sa)';
KBess  = 'at*gamma(nu+3/2)/gamma(nu)/pi^1.5./((1+sa.^2).^(nu+1.5))';
linear = 'N/A'; 
spline = 'N/A'; 
HoleS  = 'N/A'; 
HoleC  = 'N/A'; 
chris  = 'N/A'; 
wand0  = 'N/A'; 
wand1  = 'N/A'; 
wand2  = 'N/A'; 
Bohman = 'N/A'; 

Spec = char(nugget, expon, gaus, spher, cubic, penta,  cauchy, KBess,...
        linear, spline, HoleS, HoleC, chris, wand0, wand1, wand2, Bohman);

% If the model is anisotropic, consider both positive and negative frequencies
if ~iso
    ss2 = [-s(end:-1:1), s];  % Consider the negative side as well
    [s1, s2, s3] = ndgrid(s, ss2, ss2);
    s = [s1(:), s2(:), s3(:)];
    clear s1 s2 s3
end

% Initialize the spectral density matrix
fs = zeros(length(s), nvar, nvar);

% Loop through each pair (i, j) and compute the spectral densities
for i = 1:nvar
    for j = i:nvar
        for k = 1:size(model{i,j}, 1)
            % Calculate 'at' and 'sa' for isotropic and anisotropic models
            if iso
                at = model{i,j}(k, 2).^3;  % For isotropic, 'at' is a cubic power of the range
                sa = s * model{i,j}(k, 2);  % 'sa' is the frequency vector scaled by the range
            else
                [~, rot] = trans([0 0 0], model{i,j}, k);  % Rotation matrix for anisotropic model
                at = model{i,j}(k, 2) * model{i,j}(k, 3) * model{i,j}(k, 4);  % 'at' is the product of the three ranges
                sa = s * rot;  % Apply the rotation to the frequency vector
                sa = sqrt(sa.^2 * ((model{i,j}(k, 2:4)') .^ 2));  % Adjust for anisotropy
            end

            % Special treatment for numerically unstable models (penta)
            if model{i,j}(k, 1) == 6  % Unstable penta
                sa = max(sa, 0.1);  % Below 0.1, the model becomes unstable
            end
            if model{i,j}(k, 1) == 5  % Unstable penta
                sa = max(sa, 0.03);  % Below 0.03, the model becomes unstable
            end
            nu = shape{i,j}(k);
            % Evaluate the spectral density for the given model and shape
            f = eval(Spec(model{i,j}(k, 1), :));
            fs(:, i, j) = fs(:, i, j) + f * c{i,j}(k);  % Accumulate the spectral density
        end
        fs(:, j, i) = fs(:, i, j);  % Ensure symmetry
    end
end

% Permute fs for easier handling of the frequency data
fs = permute(fs, [2, 3, 1]);

% Check the eigenvalues of the spectral density matrix at each frequency
for i = 1:size(s, 1)
    ei = min(eig(fs(:,:,i)));  % Find the minimum eigenvalue at this frequency
    if ei < -1e-07 % If the minimum eigenvalue is below the threshold, the model is not admissible
        icode = 0;  % Mark the model as inadmissible
        disp(['Negative eigenvalue at frequency ', num2str(s(i)), ': ', num2str(ei)]);
        return  % Stop and return the result
    end
end

function icode = check_nugget(model, c)
% Function to check that the matrix of nugget coefficients is positive
% definite (necessary condition).
%
% Inputs:
%   model: A nvar x nvar cell array containing model parameters for each pair (i,j).
%   c: A nvar x nvar cell array containing sill parameters for each pair (i,j).
%
% Output:
%   icode: 1 if the nugget matrix is positive semi-definite, 0 otherwise.

% Number of variables (nugget matrix size)
nvar = size(model, 1);

% Initialize nugget matrix
nug = zeros(nvar, nvar);

% Iterate through the upper triangular part of the model matrix
for i = 1:nvar
    for j = i:nvar
        % Get the size of the nugget coefficient array for the current pair
        num_coeffs = size(model{i,j}, 1);

        % Check for nugget values in the model matrix
        for k = 1:num_coeffs
            if model{i,j}(k, 1) == 1  % If nugget condition is met
                nug(i,j) = c{i,j}(k);  % Set nugget coefficient
                nug(j,i) = nug(i,j);    % Ensure symmetry
            end
        end
    end
end

% Check if the nugget matrix is positive semi-definite
eig_values = eig(nug);
ei = min(eig_values);  % Get the smallest eigenvalue

% Return 1 if the matrix is positive semi-definite, otherwise 0
icode = ei >= 0;




