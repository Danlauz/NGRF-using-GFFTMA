function [dataSim] = ApplyMask(x0Ref, x0Sim, dataSim, nx, ny, nz)
    % ApplyMask - Masks the simulated data outside the convex hull of reference points
    %
    % Inputs:
    %   x0Ref  - (KxM) Matrix of reference points defining the boundary (X, Y or X, Y, Z)
    %   x0Sim  - (NxM) Matrix of simulation grid points (X, Y or X, Y, Z)
    %   dataSim - (Nx1 or NxM) Simulated data values corresponding to x0Sim
    %
    % Outputs:
    %   dataSim - Same size as input, but masked (NaN) outside the boundary

    % Find valid (non-NaN) data points
    idx = ~isnan(dataSim);

    if nz == 1 || isempty(nz)
        % 2D Case: X-Y plane
        % Compute the boundary of the valid points
        k = boundary(x0Ref(:,2), x0Ref(:,1), 0.4); % Adjust shrink factor for concave boundary
        
        % Reshape the simulation grid into a meshgrid format
        Xq = reshape(x0Sim(:,1), [nx ny]);
        Yq = reshape(x0Sim(:,2), [nx ny]);
        
        % Create a mask: Determine which points are inside the reference boundary
        [in, ~] = inpolygon(Xq, Yq, x0Ref(k,2), x0Ref(k,1));
        
        % Apply the mask: Set values outside the boundary to NaN
        dataSim(~in, :, :) = NaN;

    elseif nz > 3
        % 3D Case: X-Y-Z space
        % Compute the boundary of the valid points in 3D
        k = boundary(x0Sim(idx,1), x0Sim(idx,2), x0Sim(idx,3), 0.4); % Adjust shrink factor for concave boundary
        
        % Reshape the simulation grid into a 3D format
        Xq = reshape(x0Sim(:,1), [nx ny nz]);
        Yq = reshape(x0Sim(:,2), [nx ny nz]);
        Zq = reshape(x0Sim(:,3), [nx ny nz]);
        
        % Create a mask: Determine which points are inside the reference boundary
        [in, ~] = inpolygon(Xq, Yq, Zq, x0Ref(k,2), x0Ref(k,1), x0Ref(k,3));
        
        % Apply the mask: Set values outside the boundary to NaN
        dataSim(~in, :, :) = NaN;

    else
        % Error handling: Unexpected number of dimensions
        warning('Error in dimension: Expected 2D (X-Y) or 3D (X-Y-Z) data.');
    end
end



