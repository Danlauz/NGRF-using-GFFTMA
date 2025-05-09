function saveGridData(nx, ny, nz, sx, sy, sz, ox, oy, oz, data, filename)
    % This function saves the grid configuration and matrix data into the specified file
    %
    % Inputs:
    % - nx, ny, nz: Number of cells in the x, y, and z directions
    % - sx, sy, sz: Cell size in the x, y, and z directions
    % - ox, oy, oz: Origin coordinates in the x, y, and z directions
    % - data: A matrix of values to be saved (e.g., V1 and V2)
    % - filename: The name of the file to save the data to
    
    % Open the file for writing (create the file if it doesn't exist)
    fid = fopen(filename, 'wt');
    
    % Check if the file was opened successfully
    if fid == -1
        error('Unable to open file for writing');
    end
    
    % Write the grid parameters and cell sizes
    fprintf(fid, '# GRID - NUMBER OF CELLS\n');
    fprintf(fid, '# NX %d\n', nx);
    fprintf(fid, '# NY %d\n', ny);
    fprintf(fid, '# NZ %d\n', nz);
    
    fprintf(fid, '# GRID - CELL SIZE\n');
    fprintf(fid, '# SX %d\n', sx);
    fprintf(fid, '# SY %d\n', sy);
    fprintf(fid, '# SZ %d\n', sz);
    
    fprintf(fid, '# GRID - ORIGIN (bottom-lower-left corner)\n');
    fprintf(fid, '# OX %d\n', ox);
    fprintf(fid, '# OY %d\n', oy);
    fprintf(fid, '# OZ %d\n', oz);
    
    % Write the grid filling and sorting information
    fprintf(fid, '# GRID - FILLING\n');
    fprintf(fid, '# SORTING +X+Y+Z\n');
    
    % Write matrix values header (V1 V2)
    fprintf(fid, 'V1 V2\n');
    
    % Check if data has two columns (V1 and V2)
    [numRows, numCols] = size(data);
    if numCols == 2
        % Save the matrix data row by row
        for i = 1:numRows
            fprintf(fid, '%.6f %.6f\n', data(i, 1), data(i, 2));
        end
    else
        error('Data must have exactly two columns corresponding to V1 and V2');
    end
    
    % Close the file
    fclose(fid);
    
    disp('Data saved successfully!');
end

