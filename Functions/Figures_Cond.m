%% Visualisation of one realization (Fig1: MNGRF; Fig2: MGRF);
nvar = size(datasimNG_UC, 3); % Number of variables
% Grid
x0 = grille2(1, nx, dx, 1, ny, dy);

sample =100;

ns = 1;
% Non-Gaussian fields
FigH = figure('Position', get(0, 'Screensize'));
for i = 1:nvar

    subplot(5, nvar, i); % Arrange subplots in a single row
    imagesc(reshape(dataRefNG(:, ns, i), [nx ny])); % Reshape and display
    colormap 'jet'; % Set colormap
    colorbar(); % Add colorbar
    clim([-3.5 3.5]); % Set color limits
    cbar =colorbar();
    set(cbar, 'Ticks', -3:1:3); % Set custom ticks
    set(gca, 'YDir', 'normal', 'FontSize', 12); % Adjust axis and font size
    axis equal; % Ensure equal scaling for x and y axes
    axis tight; % Remove extra space around the plot
    xlabel('X-axis', 'FontSize', 12);
    ylabel('Y-axis', 'FontSize', 12);
    title(['V', num2str(i)], 'FontWeight', 'bold', 'FontSize', 12); % Add subplot titles

    X0 = paramA.HD{i}(:,1:2);
    LocData = (X0(:,2)-1)*nx + X0(:,1);
    subplot(5, nvar, i+nvar*1); % Arrange subplots in a single row
    imagesc(reshape(datasimNG_C_A(:, ns, i), [nx ny])); % Reshape and display
    hold on
    plot(x0(LocData,2),x0(LocData,1),'ok',MarkerSize=2)
    colormap 'jet'; % Set colormap
    colorbar(); % Add colorbar
    clim([-3.5 3.5]); % Set color limits
    cbar =colorbar();
    set(cbar, 'Ticks', -3:1:3); % Set custom ticks
    set(gca, 'YDir', 'normal', 'FontSize', 12); % Adjust axis and font size
    axis equal; % Ensure equal scaling for x and y axes
    axis tight; % Remove extra space around the plot
    xlabel('X-axis', 'FontSize', 12);
    ylabel('Y-axis', 'FontSize', 12);
    title(['V', num2str(i)], 'FontWeight', 'bold', 'FontSize', 12); % Add subplot titles
    
    X0 = paramB.HD{i}(:,1:2);
    LocData = (X0(:,2)-1)*nx + X0(:,1);
    subplot(5, nvar, i+nvar*2); % Arrange subplots in a single row
    imagesc(reshape(datasimNG_C_B(:, ns, i), [nx ny])); % Reshape and display
    hold on
    plot(x0(LocData,2),x0(LocData,1),'ok',MarkerSize=2)
    colormap 'jet'; % Set colormap
    colorbar(); % Add colorbar
    clim([-3.5 3.5]); % Set color limits
    cbar =colorbar();
    set(cbar, 'Ticks', -3:1:3); % Set custom ticks
    set(gca, 'YDir', 'normal', 'FontSize', 12); % Adjust axis and font size
    axis equal; % Ensure equal scaling for x and y axes
    axis tight; % Remove extra space around the plot
    xlabel('X-axis', 'FontSize', 12);
    ylabel('Y-axis', 'FontSize', 12);
    title(['V', num2str(i)], 'FontWeight', 'bold', 'FontSize', 12); % Add subplot titles
    
    X0 = paramC.HD{i}(:,1:2);
    LocData = (X0(:,2)-1)*nx + X0(:,1);
    subplot(5, nvar, i+nvar*3); % Arrange subplots in a single row
    imagesc(reshape(datasimNG_C_C(:, ns, i), [nx ny])); % Reshape and display
    hold on
    plot(x0(LocData,2),x0(LocData,1),'ok',MarkerSize=2)
    colormap 'jet'; % Set colormap
    colorbar(); % Add colorbar
    clim([-3.5 3.5]); % Set color limits
    cbar =colorbar();
    set(cbar, 'Ticks', -3:1:3); % Set custom ticks
    set(gca, 'YDir', 'normal', 'FontSize', 12); % Adjust axis and font size
    axis equal; % Ensure equal scaling for x and y axes
    axis tight; % Remove extra space around the plot
    xlabel('X-axis', 'FontSize', 12);
    ylabel('Y-axis', 'FontSize', 12);
    title(['V', num2str(i)], 'FontWeight', 'bold', 'FontSize', 12); % Add subplot titles
    
    X0 = paramD.HD{i}(:,1:2);
    LocData = (X0(:,2)-1)*nx + X0(:,1);
    subplot(5, nvar, i+nvar*4); % Arrange subplots in a single row
    imagesc(reshape(datasimNG_C_D(:, ns, i), [nx ny])); % Reshape and display
    hold on
    plot(x0(LocData,2),x0(LocData,1),'ok',MarkerSize=2)
    colormap 'jet'; % Set colormap
    colorbar(); % Add colorbar
    clim([-3.5 3.5]); % Set color limits
    cbar =colorbar();
    set(cbar, 'Ticks', -3:1:3); % Set custom ticks
    set(gca, 'YDir', 'normal', 'FontSize', 12); % Adjust axis and font size
    axis equal; % Ensure equal scaling for x and y axes
    axis tight; % Remove extra space around the plot
    xlabel('X-axis', 'FontSize', 12);
    ylabel('Y-axis', 'FontSize', 12);
    title(['V', num2str(i)], 'FontWeight', 'bold', 'FontSize', 12); % Add subplot titles
end

% Adjust the figure size to fit the content tightly
set(gcf, 'Position', [100, 100, 800, 1000]); % Adjust size as needed

saveas(FigH, ['Figures/Case13.png'],'png');

%% Variogram
nbsimul=min(sample, size(datasimNG_UC,2));


% Define parameters for the GeoStatFFT function
categ = 0;      % Category setting (0 for continuous)
display = 0;    % Display option
rank = 0;       % Transform data to rank distribution using ECDF
icode = 1;      % Code for variogram in GeoStatFFT

% Initialize variable to store results for each simulation
ghNG_omni_C = cell(nbsimul, 1); lagNG_omni_C = cell(nbsimul, 1);
ghNG_omni_UC = cell(nbsimul, 1); lagNG_omni_UC = cell(nbsimul, 1);
% Compute omnidirectional variogram for each realization
for i = 1:nbsimul
    [ghNG_C, nhNG_C] = GeoStatFFT(x0, datasimNG_C_B(:, i, :), icode, categ, display, rank);
    [ghNG_UC, nhNG_UC] = GeoStatFFT(x0, datasimNG_UC(:, i, :), icode, categ, display, rank);

    ang=0;
    tol_ang = 360;
    max_dist = round(2*max(nx,ny)/5); nbdist= max_dist;
    dist = [(0:nbdist-1);(1:nbdist)]'*(max_dist/nbdist);
    [ghNG_omni_C{i}, ~, lagNG_omni_C{i}] = GeoStatFFT_ndir(ghNG_C, nhNG_C, dist, ang, tol_ang);
    [ghNG_omni_UC{i}, ~, lagNG_omni_UC{i}] = GeoStatFFT_ndir(ghNG_UC, nhNG_UC, dist, ang, tol_ang);
end

% Loading Gaussian Variogram
load('BaseCaseVar.mat')

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
        meansNG_C = zeros(length(lagNG_omni_C{1}{k, j}), nbsimul);
        meansNG_UC = zeros(length(lagNG_omni_UC{1}{k, j}), nbsimul);

        % Fill quantile ranges for Gaussian fields
        for i = 1:nbsimul
            meansG(:, i) = ghG_omni{i}{k, j};
            meansNG_C(:, i) = ghNG_omni_C{i}{k, j};
            meansNG_UC(:, i) = ghNG_omni_UC{i}{k, j};
        end

        % Compute quantile ranges
        lowerBound = quantile(meansG, 0.05, 2);
        upperBound = quantile(meansG, 0.95, 2);

        % Plot the shaded confidence interval for Gaussian fields
        fill([lagG_omni{1}{k, j}; flipud(lagG_omni{1}{k, j})], ...
            [lowerBound; flipud(upperBound)], [0 0.8 0], ...
            'EdgeColor', 'none', 'FaceAlpha', 0.3, 'DisplayName', 'C.I. 95% Gaussian Fields');


        % Compute quantile ranges
        lowerBound = quantile(meansNG_C, 0.05, 2);
        upperBound = quantile(meansNG_C, 0.95, 2);

        % Plot the shaded confidence interval for Non- Gaussian fields
        fill([lagNG_omni_C{1}{k, j}; flipud(lagNG_omni_C{1}{k, j})], ...
            [lowerBound; flipud(upperBound)], [0 0 0.9], ...
            'EdgeColor', 'none', 'FaceAlpha', 0.3, 'DisplayName', 'C.I. 95% non-Gaussian Fields (C)');

        % Compute quantile ranges
        lowerBound = quantile(meansNG_UC, 0.05, 2);
        upperBound = quantile(meansNG_UC, 0.95, 2);

        % Plot the shaded confidence interval for Non- Gaussian fields
        fill([lagNG_omni_UC{1}{k, j}; flipud(lagNG_omni_UC{1}{k, j})], ...
            [lowerBound; flipud(upperBound)], [0.9 0 0], ...
            'EdgeColor', 'none', 'FaceAlpha', 0.3, 'DisplayName', 'C.I. 95% non-Gaussian Fields (UC)');


        % Plot the average for Gaussian fields
        plot(lagG_omni{1}{k, j}, mean(meansG, 2), '-', ...
            'LineWidth', 6, 'Color', [0 0.5 0], 'DisplayName', 'Average Across Gaussian Realizations');

        % Plot the average for non-Gaussian fields
        plot(lagNG_omni_C{1}{k, j}, mean(meansNG_C, 2), '-', ...
            'LineWidth', 6, 'Color', [0 0 0.6], 'DisplayName', 'Average Across non-Gaussian Realizations');

        % Plot the average for non-Gaussian fields
        plot(lagNG_omni_UC{1}{k, j}, mean(meansNG_UC, 2), '-', ...
            'LineWidth', 6, 'Color', [0.6 0 0], 'DisplayName', 'Average Across non-Gaussian Realizations');


        % Compute and plot the theoretical covariance function
        if length(model{k, j})==2
            K = c{k, j} - covar((0:nx)', 0, model{k, j}, c{k, j}, nu{k, j});
        else
            K = c{k, j} - covar([(0:nx)' zeros(nx+1,1) zeros(nx+1,1)] , [0 0 0], model{k, j}, c{k, j}, nu{k, j});
        end
        plot(0:nx, K, '--', 'LineWidth', 6, 'Color', [0 0 0], 'DisplayName', 'Theoretical');

        % Labeling and formatting
        title(['V', num2str(k), ' vs. V', num2str(j)], 'FontWeight', 'bold', 'FontSize', 28);
        xlim([0, round(2 * nx / 5)]);
        ylim([0, 1.5]);
        xlabel('Distance', 'FontSize', 28);
        ylabel('Var. Value', 'FontSize', 28);
        set(gca, 'FontSize', 28);
        %legend('show');
        grid on;
        hold off;

    end
end

saveas(FigH, ['Figures/Case' num2str(cas) '_Var.png'],'png');
%% Rank Asymmetry
nbsimul=min(sample, size(datasimNG_C_B,2));

% Grid
x0 = grille2(1, nx, dx, 1, ny, dy);

% Define parameters for the GeoStatFFT function
categ = 0;      % Category setting (0 for continuous)
display = 0;    % Display option
rank = 1;       % Transform data to rank distribution using ECDF
icode = 9;      % Code for variogram in GeoStatFFT

% Initialize variable to store results for each simulation
ghNG_omni_C = cell(nbsimul, 1); lagNG_omni_C = cell(nbsimul, 1);
ghNG_omni_UC = cell(nbsimul, 1); lagNG_omni_UC = cell(nbsimul, 1);
% Compute omnidirectional variogram for each realization
for i = 1:nbsimul
    [ghNG_C, nhNG_C] = GeoStatFFT(x0, datasimNG_C_B(:, i, :), icode, categ, display, rank);
    [ghNG_UC, nhNG_UC] = GeoStatFFT(x0, datasimNG_UC(:, i, :), icode, categ, display, rank);

    ang=0;
    tol_ang = 360;
    max_dist = round(2*max(nx,ny)/5); nbdist= max_dist;
    dist = [(0:nbdist-1);(1:nbdist)]'*(max_dist/nbdist);
    [ghNG_omni_C{i}, ~, lagNG_omni_C{i}] = GeoStatFFT_ndir(ghNG_C, nhNG_C, dist, ang, tol_ang);
    [ghNG_omni_UC{i}, ~, lagNG_omni_UC{i}] = GeoStatFFT_ndir(ghNG_UC, nhNG_UC, dist, ang, tol_ang);
end

% Loading Gaussian Rank Asymetry
load('BaseCaseRankAss.mat')

% Plot Rank Asy. for each pair of variables
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
        meansNG_C = zeros(length(lagNG_omni_C{1}{k, j}), nbsimul);
        meansNG_UC = zeros(length(lagNG_omni_UC{1}{k, j}), nbsimul);

        % Fill quantile ranges for Gaussian fields
        for i = 1:nbsimul
            meansG(:, i) = ghG_omni{i}{k, j};
            meansNG_C(:, i) = ghNG_omni_C{i}{k, j};
            meansNG_UC(:, i) = ghNG_omni_UC{i}{k, j};
        end

        % Compute quantile ranges
        lowerBound = quantile(meansG, 0.05, 2);
        upperBound = quantile(meansG, 0.95, 2);

        % Plot the shaded confidence interval for Gaussian fields
        fill([lagG_omni{1}{k, j}; flipud(lagG_omni{1}{k, j})], ...
            [lowerBound; flipud(upperBound)], [0 0.8 0], ...
            'EdgeColor', 'none', 'FaceAlpha', 0.3, 'DisplayName', 'C.I. 95% Gaussian Fields');


        % Compute quantile ranges
        lowerBound = quantile(meansNG_C, 0.05, 2);
        upperBound = quantile(meansNG_C, 0.95, 2);

        % Plot the shaded confidence interval for Non- Gaussian fields
        fill([lagNG_omni_C{1}{k, j}; flipud(lagNG_omni_C{1}{k, j})], ...
            [lowerBound; flipud(upperBound)], [0 0 0.9], ...
            'EdgeColor', 'none', 'FaceAlpha', 0.3, 'DisplayName', 'C.I. 95% non-Gaussian Fields (C)');

        % Compute quantile ranges
        lowerBound = quantile(meansNG_UC, 0.05, 2);
        upperBound = quantile(meansNG_UC, 0.95, 2);

        % Plot the shaded confidence interval for Non- Gaussian fields
        fill([lagNG_omni_UC{1}{k, j}; flipud(lagNG_omni_UC{1}{k, j})], ...
            [lowerBound; flipud(upperBound)], [0.9 0 0], ...
            'EdgeColor', 'none', 'FaceAlpha', 0.3, 'DisplayName', 'C.I. 95% non-Gaussian Fields (UC)');


        % Plot the average for Gaussian fields
        plot(lagG_omni{1}{k, j}, mean(meansG, 2), '-', ...
            'LineWidth', 6, 'Color', [0 0.5 0], 'DisplayName', 'Average Across Gaussian Realizations');

        % Plot the average for non-Gaussian fields
        plot(lagNG_omni_C{1}{k, j}, mean(meansNG_C, 2), '-', ...
            'LineWidth', 6, 'Color', [0 0 0.6], 'DisplayName', 'Average Across non-Gaussian Realizations');

        % Plot the average for non-Gaussian fields
        plot(lagNG_omni_UC{1}{k, j}, mean(meansNG_UC, 2), '-', ...
            'LineWidth', 6, 'Color', [0.6 0 0], 'DisplayName', 'Average Across non-Gaussian Realizations');

        % Compute and plot the theoretical covariance function
        plot([0 nx], [0 0], '--', 'LineWidth', 6, 'Color', [0 0 0], 'DisplayName', 'Theoretical');


        % Labeling and formatting
        title(['V', num2str(k), ' vs. V', num2str(j)], 'FontWeight', 'bold', 'FontSize', 28);
        xlim([0, round(2 * nx / 5)]);
        ylim([-0.05, 0.05]);
        xlabel('Distance', 'FontSize', 28);
        ylabel('Var. Value', 'FontSize', 28);
        set(gca, 'FontSize', 28);
        %legend('show');
        grid on;
        hold off;

    end
end
saveas(FigH, ['Figures/Case' num2str(cas) '_RankAss.png'],'png');

%% Directional Asymmetry
nbsimul=min(sample, size(datasimNG_C_B,2));

% Grid
x0 = grille2(1, nx, dx, 1, ny, dy);

% Define parameters for the GeoStatFFT function
categ = 0;      % Category setting (0 for continuous)
display = 0;    % Display option
rank = 1;       % Transform data to rank distribution using ECDF
icode = 8;      % Code for variogram in GeoStatFFT

% Initialize variable to store results for each simulation
ghG = cell(nbsimul, 1); ghNG_C = cell(nbsimul, 1); ghNG_UC = cell(nbsimul, 1);
% Compute omnidirectional variogram for each realization
for i = 1:nbsimul
    ghNG_C{i} = GeoStatFFT(x0, datasimNG_C_B(:, i, :), icode, categ, display, rank);
    ghNG_UC{i} = GeoStatFFT(x0, datasimNG_UC(:, i, :), icode, categ, display, rank);
end

% Loading Gaussian Directional Asymetry
load('BaseCaseDirAss.mat')

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
        meansG = zeros(length(ghG_dirX{1}{k, j}), nbsimul);
        meansNG_C = zeros(length(ghNG_C{1}{k, j}(nx:end, ny)), nbsimul);
        meansNG_UC = zeros(length(ghNG_UC{1}{k, j}(nx:end, ny)), nbsimul);

        % Fill quantile ranges for Gaussian fields
        for i = 1:nbsimul
            meansG(:, i) = ghG_dirX{i}{k, j};
            meansNG_C(:, i) = ghNG_C{i}{k, j}(nx, ny:end);
            meansNG_UC(:, i) = ghNG_UC{i}{k, j}(nx, ny:end);
        end

        % Compute quantile ranges
        lowerBound = quantile(meansG, 0.05, 2);
        upperBound = quantile(meansG, 0.95, 2);

        % Plot the shaded confidence interval for Gaussian fields
        x = (0:nx-1)';
        fill([x; flipud(x)], ...
            [lowerBound; flipud(upperBound)], [0 0.9 0], ...
            'EdgeColor', 'none', 'FaceAlpha', 0.3, 'DisplayName', 'C.I. 95% Gaussian Fields');

        % Compute quantile ranges
        lowerBound = quantile(meansNG_C, 0.05, 2);
        upperBound = quantile(meansNG_C, 0.95, 2);

        % Plot the shaded confidence interval for Gaussian fields
        x = (0:nx-1)';
        fill([x; flipud(x)], ...
            [lowerBound; flipud(upperBound)], [0 0 0.9], ...
            'EdgeColor', 'none', 'FaceAlpha', 0.3, 'DisplayName', 'C.I. 95% non-Gaussian Fields');

        % Compute quantile ranges
        lowerBound = quantile(meansNG_UC, 0.05, 2);
        upperBound = quantile(meansNG_UC, 0.95, 2);

        % Plot the shaded confidence interval for Gaussian fields
        x = (0:nx-1)';
        fill([x; flipud(x)], ...
            [lowerBound; flipud(upperBound)], [0.9 0 0], ...
            'EdgeColor', 'none', 'FaceAlpha', 0.3, 'DisplayName', 'C.I. 95% non-Gaussian Fields');


        % Plot the average variogram for Gaussian fields
        plot(x, mean(meansG, 2), '-', ...
            'LineWidth', 6, 'Color', [0 0.8 0], 'DisplayName', 'Average Across Realizations');

        % Plot the average variogram for Gaussian fields
        plot(x, mean(meansNG_C, 2), '-', ...
            'LineWidth', 6, 'Color', [0 0 0.8], 'DisplayName', 'Average Across Realizations');

        % Plot the average variogram for Gaussian fields
        plot(x, mean(meansNG_UC, 2), '-', ...
            'LineWidth', 6, 'Color', [0.8 0 0], 'DisplayName', 'Average Across Realizations');


        % Compute and plot the theoretical covariance function
        plot([0 nx], [0 0], '--', 'LineWidth', 6, 'Color', [0 0 0], 'DisplayName', 'Theoretical');

        % Labeling and formatting
        title(['V', num2str(k), ' vs. V', num2str(j)], 'FontWeight', 'bold', 'FontSize', 28);
        xlim([0, round(1 * nx / 5)]);
        ylim([-0.05, 0.05]);
        xlabel('Distance', 'FontSize', 28);
        ylabel('Dir. Asy. Value', 'FontSize', 28);
        set(gca, 'FontSize', 28);
        %legend('show');
        grid on;
        hold off;
    end
end
saveas(FigH, ['Figures/Case' num2str(cas) '_DirAss.png'],'png');