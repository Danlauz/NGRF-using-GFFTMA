%% Visualisation of one realization (Fig1: MGRF);
FigH = figure('Position', get(0, 'Screensize'));
nvar = size(datasimG, 3); % Number of variables
for i = 1:nvar
    subplot(1, nvar, i); % Arrange subplots in a single row
    imagesc(reshape(datasimG(:, 1, i), [nx ny])); % Reshape and display
    colormap 'jet'; % Set colormap
    colorbar(); % Add colorbar
    clim([-3.5 3.5]); % Set color limits
    cbar =colorbar();
    set(cbar, 'Ticks', -3:1:3); % Set custom ticks
    set(gca, 'YDir', 'normal', 'FontSize', 32); % Adjust axis and font size
    axis equal; % Ensure equal scaling for x and y axes
    axis tight; % Remove extra space around the plot
    xlabel('X-axis', 'FontSize', 32);
    ylabel('Y-axis', 'FontSize', 32);
    title(['V', num2str(i)], 'FontWeight', 'bold', 'FontSize', 32); % Add subplot titles
end

% Adjust the figure size to fit the content tightly
set(gcf, 'Position', [100, 100, 1600, 625]); % Adjust size as needed

save('BaseCaseSim.mat','datasimG')
saveas(FigH, ['Figures/Case' num2str(cas) '_Real.png'],'png');
sample = 100;
%% Variogram
nbsimul=min(sample, size(datasimG,2));

% Grid
x0 = grille2(1, nx, dx, 1, ny, dy);

% Define parameters for the GeoStatFFT function
categ = 0;      % Category setting (0 for continuous)
display = 0;    % Display option
rank = 0;       % Transform data to rank distribution using ECDF
icode = 1;      % Code for variogram in GeoStatFFT

% Initialize variable to store results for each simulation
ghG = cell(nbsimul, 1);

% Compute omnidirectional variogram for each realization
for i = 1:nbsimul
    [ghG, nhG] = GeoStatFFT(x0, datasimG(:, i, :), icode, categ, display, rank);

    ang=0;
    tol_ang = 360;
    max_dist = round(2*max(nx,ny)/5); nbdist= max_dist;
    dist = [(0:nbdist-1);(1:nbdist)]'*(max_dist/nbdist);
    [ghG_omni{i}, ~, lagG_omni{i}] = GeoStatFFT_ndir(ghG, nhG, dist, ang, tol_ang);
end

% Plot variogram for each pair of variables
FigH = figure('Position', get(0, 'Screensize'));
nvar = size(model, 1); % Number of variables
a = 0; % Counter for subplot indexing

% Iterate over variable pairs
for k = 1:nvar
    for j = k:nvar
        a = a + 1;
        subplot(nvar, nvar, (k-1)*nvar + j); % Adjust subplot grid for clarity
        hold on;

        % Initialize matrices for quantile computation
        meansG = zeros(length(lagG_omni{1}{k, j}), nbsimul);

        % Fill quantile ranges for Gaussian fields
        for i = 1:nbsimul
            meansG(:, i) = ghG_omni{i}{k, j};
        end

        % Compute quantile ranges
        lowerBound = quantile(meansG, 0.05, 2);
        upperBound = quantile(meansG, 0.95, 2);

        % Plot the shaded confidence interval for Gaussian fields
        fill([lagG_omni{1}{k, j}; flipud(lagG_omni{1}{k, j})], ...
            [lowerBound; flipud(upperBound)], [0 0.8 0], ...
            'EdgeColor', 'none', 'FaceAlpha', 0.3, 'DisplayName', 'C.I. 95% Gaussian Fields');

        % Plot the average variogram for Gaussian fields
        plot(lagG_omni{1}{k, j}, mean(meansG, 2), '-', ...
            'LineWidth', 6, 'Color', [0 0.5 0], 'DisplayName', 'Average Across Realizations');

        % Compute and plot the theoretical covariance function
        K = c{k, j} - covar((0:nx)', 0, model{k, j}, c{k, j}, nu{k, j});
        plot(0:nx, K, '--', 'LineWidth', 6, 'Color', [0 0 0], 'DisplayName', 'Theoretical');

        % Labeling and formatting
        title(['V', num2str(k), ' vs. V', num2str(j)], 'FontWeight', 'bold', 'FontSize', 28);
        xlim([0, round(2 * nx / 5)]);
        ylim([0, 1.5]);
        xlabel('Distance', 'FontSize', 28);
        ylabel('Var. Value', 'FontSize', 28);
        set(gca, 'FontSize', 28);
        grid on;
        hold off;
    end
end
set(gcf, 'Position', [100, 100, 1600, 625]); % Adjust size as needed

save('BaseCaseVar.mat','ghG_omni', 'lagG_omni')
saveas(FigH, ['Figures/Case' num2str(cas) '_Var.png'],'png');
%% Rank Asymmetry
nbsimul=min(sample, size(datasimG,2));

% Grid
x0 = grille2(1, nx, dx, 1, ny, dy);

% Define parameters for the GeoStatFFT function
categ = 0;      % Category setting (0 for continuous)
display = 0;    % Display option
rank = 1;       % Transform data to rank distribution using ECDF
icode = 9;      % Code for variogram in GeoStatFFT

% Initialize variable to store results for each simulation
ghG = cell(nbsimul, 1);

% Compute omnidirectional variogram for each realization
for i = 1:nbsimul
    [ghG, nhG] = GeoStatFFT(x0, datasimG(:, i, :), icode, categ, display, rank);

    ang=0;
    tol_ang = 360;
    max_dist = round(2*max(nx,ny)/5); nbdist= max_dist;
    dist = [(0:nbdist-1);(1:nbdist)]'*(max_dist/nbdist);
    [ghG_omni{i}, ~, lagG_omni{i}] = GeoStatFFT_ndir(ghG, nhG, dist, ang, tol_ang);
end

% Plot variogram for each pair of variables
FigH = figure('Position', get(0, 'Screensize'));
nvar = size(model, 1); % Number of variables
a = 0; % Counter for subplot indexing

% Iterate over variable pairs
for k = 1:nvar
    for j = k:nvar
        a = a + 1;
        subplot(nvar, nvar, (k-1)*nvar + j); % Adjust subplot grid for clarity
        hold on;

        % Initialize matrices for quantile computation
        meansG = zeros(length(lagG_omni{1}{k, j}), nbsimul);

        % Fill quantile ranges for Gaussian fields
        for i = 1:nbsimul
            meansG(:, i) = ghG_omni{i}{k, j};
        end

        % Compute quantile ranges
        lowerBound = quantile(meansG, 0.05, 2);
        upperBound = quantile(meansG, 0.95, 2);

        % Plot the shaded confidence interval for Gaussian fields
        fill([lagG_omni{1}{k, j}; flipud(lagG_omni{1}{k, j})], ...
            [lowerBound; flipud(upperBound)], [0 0.8 0], ...
            'EdgeColor', 'none', 'FaceAlpha', 0.3, 'DisplayName', 'C.I. 95% Gaussian Fields');

        % Plot the average variogram for Gaussian fields
        plot(lagG_omni{1}{k, j}, mean(meansG, 2), '-', ...
            'LineWidth', 6, 'Color', [0 0.5 0], 'DisplayName', 'Average Across Realizations');

        % Compute and plot the theoretical covariance function
        plot([0 nx], [0 0], '--', 'LineWidth', 6, 'Color', [0 0 0], 'DisplayName', 'Theoretical');

        % Labeling and formatting
        title(['V', num2str(k), ' vs. V', num2str(j)], 'FontWeight', 'bold', 'FontSize', 28);
        xlim([0, round(1 * nx / 5)]);
        ylim([-0.05 0.05]);
        xlabel('Distance', 'FontSize', 28);
        ylabel('Rank Asy. Value', 'FontSize', 28);
        set(gca, 'FontSize', 28);
        grid on;
        hold off;
    end
end

save('BaseCaseRankAss.mat','ghG_omni', 'lagG_omni')
saveas(FigH, ['Figures/Case' num2str(cas) '_RankAss.png'],'png');
%% Directional Asymmetry
nbsimul=min(sample, size(datasimG,2));

% Grid
x0 = grille2(1, nx, dx, 1, ny, dy);

% Define parameters for the GeoStatFFT function
categ = 0;      % Category setting (0 for continuous)
display = 0;    % Display option
rank = 1;       % Transform data to rank distribution using ECDF
icode = 8;      % Code for variogram in GeoStatFFT

% Initialize variable to store results for each simulation
ghG = cell(nbsimul, 1);

% Compute omnidirectional variogram for each realization
for i = 1:nbsimul
    ghG{i} = GeoStatFFT(x0, datasimG(:, i, :), icode, categ, display, rank);
end

% Plot variogram for each pair of variables
FigH = figure('Position', get(0, 'Screensize'));
nvar = size(model, 1); % Number of variables
a = 0; % Counter for subplot indexing

% Initialize variable to store results for each simulation
ghG_dirX = cell(nbsimul, 1);
for i = 1:nbsimul
    for k = 1:nvar
        for j = k:nvar
            ghG_dirX{i}{k,j} = ghG{i}{k, j}(nx:end, ny);
        end
    end
end
% Iterate over variable pairs
for k = 1:nvar
    for j = k:nvar
        a = a + 1;
        subplot(nvar, nvar, (k-1)*nvar + j); % Adjust subplot grid for clarity
        hold on;

        % Initialize matrices for quantile computation
        meansG = zeros(length(ghG_dirX{i}{k,j}), nbsimul);

        % Fill quantile ranges for Gaussian fields
        for i = 1:nbsimul
            meansG(:, i) = ghG_dirX{i}{k,j};
        end

        % Compute quantile ranges
        lowerBound = quantile(meansG, 0.05, 2);
        upperBound = quantile(meansG, 0.95, 2);

        % Plot the shaded confidence interval for Gaussian fields
        x = (0:nx-1)';
        fill([x; flipud(x)], ...
            [lowerBound; flipud(upperBound)], [0 0.8 0], ...
            'EdgeColor', 'none', 'FaceAlpha', 0.3, 'DisplayName', 'C.I. 95% Gaussian Fields');

        % Plot the average variogram for Gaussian fields
        plot(x, mean(meansG, 2), '-', ...
            'LineWidth', 6, 'Color', [0 0.5 0], 'DisplayName', 'Average Across Realizations');

        % Compute and plot the theoretical covariance function
        plot([0 nx], [0 0], '--', 'LineWidth', 6, 'Color', [0 0 0], 'DisplayName', 'Theoretical');

        % Labeling and formatting
        title(['V', num2str(k), ' vs. V', num2str(j)], 'FontWeight', 'bold', 'FontSize', 28);
        xlim([0, round(1 * nx / 5)]);
        ylim([-0.05, 0.05]);
        xlabel('Distance', 'FontSize', 28);
        ylabel('Dir. Asy. Value', 'FontSize', 28);
        set(gca, 'FontSize', 28);
        grid on;
        hold off;
    end
end
save('BaseCaseDirAss.mat','ghG_dirX')
saveas(FigH, ['Figures/Case' num2str(cas) '_DirAss.png'],'png');