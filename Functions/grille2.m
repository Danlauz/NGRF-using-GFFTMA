function gril = grille2(xmin, xmax, dx, ymin, ymax, dy)
    % Function to create a regular 2D grid
    % The grid spans from xmin to xmax and ymin to ymax with steps dx and dy.
    
    % Generate x and y coordinates
    x = (xmin:dx:xmax)';  % Column vector for x
    y = (ymin:dy:ymax)';  % Column vector for y
    
    % Get the number of points in each direction
    nx = length(x);
    ny = length(y);
    
    % Create the grid by combining x and y coordinates
    gril = [kron(ones(ny, 1), x), kron(y, ones(nx, 1))];
end