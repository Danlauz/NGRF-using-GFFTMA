%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 41);

% Specify sheet and range
opts.Sheet = "2mm";
opts.DataRange = "A3:AO113";

% Specify column names and types
opts.VariableNames = ["Var1", "Var2", "Var3", "Easting", "Northing", "Var6", "Elevation", "Agppm", "Al", "Asppm", "Var11", "Var12", "Var13", "Var14", "Ca", "Cdppm", "Var17", "Var18", "Cuppm", "Var20", "Var21", "Var22", "K", "Var24", "Mg", "Var26", "Var27", "Var28", "Nippm", "Var30", "Pbppm", "Var32", "Var33", "Var34", "Var35", "Var36", "Var37", "Var38", "Vppm", "Var40", "Znppm"];
opts.SelectedVariableNames = ["Easting", "Northing", "Agppm", "Al", "Asppm", "Ca", "Cdppm", "Cuppm", "K", "Mg", "Nippm", "Pbppm", "Vppm", "Znppm"];
opts.VariableTypes = ["char", "char", "char", "double", "double", "char", "double", "double", "double", "double", "char", "char", "char", "char", "double", "double", "char", "char", "double", "char", "char", "char", "double", "char", "double", "char", "char", "char", "double", "char", "double", "char", "char", "char", "char", "char", "char", "char", "double", "char", "double"];

% Specify variable properties
opts = setvaropts(opts, ["Var1", "Var2", "Var3", "Var6", "Var11", "Var12", "Var13", "Var14", "Var17", "Var18", "Var20", "Var21", "Var22", "Var24", "Var26", "Var27", "Var28", "Var30", "Var32", "Var33", "Var34", "Var35", "Var36", "Var37", "Var38", "Var40"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var1", "Var2", "Var3", "Var6", "Var11", "Var12", "Var13", "Var14", "Var17", "Var18", "Var20", "Var21", "Var22", "Var24", "Var26", "Var27", "Var28", "Var30", "Var32", "Var33", "Var34", "Var35", "Var36", "Var37", "Var38", "Var40"], "EmptyFieldRule", "auto");

% Import the data
Data = readtable("C:\Users\dany-\OneDrive - polymtl\08_Projet_EnCours\Article12_NG-FFTMA\MPS_Conditioning\Horne_OrganicSoils_ICP_AES_ARD_2mm_C11.xls", opts, "UseExcel", false);
Data = Data([1:98 108:end],:);

[x, y] = utm2ll(Data.Easting, Data.Northing, 17);
Data.Easting = x;
Data.Northing = y;

clear opts

% Convert to Lat/Lon
CoorHorne = [48.2533 -79.0136];
%% Grid resolution
data = table2array(Data);
data = data(:,[1 2 5 7]);
elementTitles = { "log(As) (ppm)",  "log(Cd) (ppm)",};

res = 0.01;
Data.Asppm = log(data(:,3));
[dataAs, gridAs, nx, ny] = grid_data(Data , 'Easting', 'Northing', 'Asppm', res);

Data.Cdppm = log(data(:,4));
[dataCd, gridCd, ~, ~] = grid_data(Data , 'Easting', 'Northing', 'Cdppm', res);

data = [dataAs.X, dataAs.Y,  dataAs.Z, dataCd.Z];

x0 = grille2(1,nx,1,1,ny,1);
%% Gaussian transformation
for i=1:2
    idx = ~isnan(data(:,2+i));
    datamean(i) = mean(data(idx ,2+i));
    datastd(i) =std(data(idx ,2+i));
    data(idx ,2+i) = (data(idx ,2+i)-datamean(i))/datastd(i);
end
%% Omnidirectional variogram (fitting was done manually)

model{1,1} = [2 44/3];   model{1,2} = [2 38/3];       
model{2,1} = model{1,2};   model{2,2} = [2 35/3]; 
 
c{1,1} = 1;        c{1,2} = 0.9;      
c{2,1} = c{1,2};   c{2,2} = 1;

nu{1,1} = 1;       nu{1,2} = 1;    
nu{2,1} = nu{1,2}; nu{2,2} = 1;  

s=[0.0005:0.0001:0.01 0.02:0.001:0.1 0.2:0.01:7 7.2:0.1:100 101:1:1000]'; % set of frenquencies to check
icode = tasc3d(model, c, nu, s);
disp(['Admissible ', num2str(mean(icode) * 100), '%']);

%% Rank Asymmetry (fitting was done manually)
nbfam = 21;

% For Non-Gaussian modelling
param.rangeX = cell(2, 2); param.nu = cell(2, 2);
param.rangeY = cell(2, 2); param.c = cell(2, 2);
param.Gaus = [false false]; param.Below = [1 1];
paramG = param;

param.rangeX{1,1} = [1 88*4]/3;          param.rangeX{1,2} = [1 76*4]/3;         
param.rangeX{2,1} = param.rangeX{1,2};  param.rangeX{2,2} = [1 70*4]/3; 

param.rangeY = cell(2, 2); 
param.rotX = cell(2, 2);
param.dir = cell(2, 1);

% For Gaussian modelling
paramG.Gaus = [true true];
paramG.rangeX{1,1} = [44 44]/3;          paramG.rangeX{1,2} = [38 38]/3;
paramG.rangeX{2,1} = param.rangeX{1,2};  paramG.rangeX{2,2} = [35 35]/3;


s=[0.0005:0.0001:0.01 0.02:0.001:0.1 0.2:0.01:7 7.2:0.1:100 101:1:1000]'; % set of frenquencies to check
icode = NG_tasc3d(model, c, nu, s, param, nbfam);
disp(['Family is admissible ', num2str(mean(icode) * 100), '% of the time.']);

%% Conditioning Data
a = 0;
for i=[1 2]
    a= a+1;
    idx = ~isnan(data(:,2+i));
    param.HD{a} = [x0(idx, 1:2) data(idx, 2+i)] ;
    paramG.HD{a} = [x0(idx, 1:2) data(idx, 2+i)] ;
end

%% Simulation
nbsimul = 100;
seed = 15451;

% Multivariate Non-Gaussian Random Fields
tic
[datasimNG, ErrNG] = GFFTMA_NG(model, c, nu, param, seed, nbsimul, nbfam, nx, 1, ny, 1);
tNG = toc;

% Multivariate Gaussian Random Fields
tic
datasimG= GFFTMA_NG(model, c, nu, paramG, seed, nbsimul, 1, nx, 1, ny, 1);
tG = toc;
Data = table2array(Data);

%% Figures
%% Tranform to log-space and original-space
for i=[1 2]
    datasimNG(:,:,i) = datasimNG(:,:,i)*datastd(i) + datamean(i);
    datasimG(:,:,i) = datasimG(:,:,i)*datastd(i) + datamean(i);
    datasimNG_Org(:,:,i) = exp(datasimNG(:,:,i));
    datasimG_Org(:,:,i) = exp(datasimG(:,:,i));
end

% Data transformation into km units
[data(:,1),data(:,2),~]=ll2utm(data(:,1), data(:,2));
[Data(:,1),Data(:,2),~]=ll2utm(Data(:,1), Data(:,2));
[CoorHorne(:,1),CoorHorne(:,2),f]=ll2utm(CoorHorne(:,1), CoorHorne(:,2));

data(:,1) = (data(:,1) -CoorHorne(:,1))/1000; data(:,2) = (data(:,2) -CoorHorne(:,2))/1000;
Data(:,1) = (Data(:,1) -CoorHorne(:,1))/1000; Data(:,2) = (Data(:,2) -CoorHorne(:,2))/1000;
CoorHorne(:,1) =0; CoorHorne(:,2)=0;

minY = min(data(:,1))-5; maxY=max(data(:,1))+5;
minX = min(data(:,2))-5; maxX=max(data(:,2))+5;

%% Figure Background
idx = ~isnan(Data(:,5)) | ~isnan(Data(:,7)) ;
Data = Data(idx,[1 2 5 7]);

figure(1)
tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Loop through the tiles

a =0;
for i = [1 2] 
    a = a +1;
    nexttile;
    Z = Data(:, i + 2); % Extract data for the current element

    % Filter out NaN values
    idx = ~isnan(Z);
    Z = Z(idx);
    % Determine color limits dynamically to account for negative values
    climRange = [min(Z) max(Z)]; 

    % Scatter plot with adjusted color mapping
    scatter(Data(idx, 2), Data(idx, 1), 40, Z, 'filled');
    hold on
    plot(CoorHorne(2), CoorHorne(1), 'pk', 'MarkerSize', 10, 'MarkerFaceColor', 'k')

    % Configure colorbar and colormap
    colorbar();
    colormap('jet'); % Change to 'balance' or 'redblue' if you prefer a zero-centered colormap
    clim(climRange); % Adjust color limits dynamically

    % Set axis limits and labels
    axis equal
    xlim([-0.1+minX, maxX + 0.1]);
    ylim([-0.1+minY, maxY + 0.1]);
    xlabel('Easting (km)');
    ylabel('Northing (km)');
    set(gca, 'FontSize', 14, 'FontUnits', 'points');
    set(gca, 'XMinorTick', 'on');
    set(gca, 'YMinorTick', 'on');

    % Add the title
    title(elementTitles{i});
end
%% Visualisation of one realization (Fig100: MNGRF; Fig100: MGRF);

% Figure of one realization
a = 0;
figure(97)
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
elementTitles = { "log(As) (ppm)",  "log(Cd) (ppm)",};
for i = [1 2]
    a= a+1;
    nexttile;

    % Apply a mask to data to eliminate extrapolation
    Z = ApplyMask(Data(:,[2 1]), data(:,[1 2]), datasimNG(:, 2, a), nx ,ny, 1);

    % Determine color limits dynamically to account for negative values
    if a ==1
        climRange = [0 5];
    else
        climRange = [-2 4];
    end   

    % Scatter plot with adjusted color mapping
    scatter(data(:, 2), data(:, 1), 20, Z, 'filled',"square");
    hold on
    plot(CoorHorne(2), CoorHorne(1), 'pk', 'MarkerSize', 10, 'MarkerFaceColor', 'k')
    hold on
    scatter(Data(:, 2), Data(:, 1), 20, 'k', 'filled', 'c');

    % Configure colorbar and colormap
    colorbar();
    colormap('jet'); % Change to 'balance' or 'redblue' if you prefer a zero-centered colormap
    clim(climRange); % Adjust color limits dynamically

    % Set axis limits and labels
    axis equal
    xlim([minX, maxX]);
    ylim([minY, maxY]);
    xlabel('Easting (km)');
    ylabel('Northing (km)');
    set(gca, 'FontSize', 14, 'FontUnits', 'points');
    set(gca, 'XMinorTick', 'on');
    set(gca, 'YMinorTick', 'on');

    % Add the title
    title(elementTitles{i});

    nexttile;

    % Apply a mask to data to eliminate extrapolation
    Z = ApplyMask(Data(:,[2 1]), data(:,[1 2]), datasimG(:, 2, a), nx ,ny, 1);

    % Determine color limits dynamically to account for negative values
    if a ==1
        climRange = [0 5];
    else
        climRange = [-2 4];
    end   

    % Scatter plot with adjusted color mapping
    scatter(data(:, 2), data(:, 1), 20, Z, 'filled',"square");
    hold on
    plot(CoorHorne(2), CoorHorne(1), 'pk', 'MarkerSize', 10, 'MarkerFaceColor', 'k')
    hold on
    scatter(Data(:, 2), Data(:, 1), 20, 'k', 'filled', 'c');

    % Configure colorbar and colormap
    colorbar();
    colormap('jet'); % Change to 'balance' or 'redblue' if you prefer a zero-centered colormap
    clim(climRange); % Adjust color limits dynamically

    % Set axis limits and labels
    axis equal
    xlim([minX, maxX]);
    ylim([minY, maxY]);
    xlabel('Easting (km)');
    ylabel('Northing (km)');
    set(gca, 'FontSize', 14, 'FontUnits', 'points');
    set(gca, 'XMinorTick', 'on');
    set(gca, 'YMinorTick', 'on');

    % Add the title
    title(elementTitles{i});
end

%% Figure of mean
a = 0;
figure(98)
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
elementTitles = { "log(As) (ppm)",  "log(Cd) (ppm)",};
for i = [1 2]
    a= a+1;
    nexttile;

    % Apply a mask to data to eliminate extrapolation
    Z = ApplyMask(Data(:,[2 1]), data(:,[1 2]), mean(datasimNG(:, :, a),2), nx ,ny, 1);

    % Determine color limits dynamically to account for negative values
    if a ==1
        climRange = [1 4];
    else
        climRange = [0 3];
    end   

    % Scatter plot with adjusted color mapping
    scatter(data(:, 2), data(:, 1), 20, Z, 'filled',"square");
    hold on
    plot(CoorHorne(2), CoorHorne(1), 'pk', 'MarkerSize', 10, 'MarkerFaceColor', 'k')
    hold on
    scatter(Data(:, 2), Data(:, 1), 20, 'k', 'filled', 'c');

    % Configure colorbar and colormap
    colorbar();
    colormap('jet'); % Change to 'balance' or 'redblue' if you prefer a zero-centered colormap
    clim(climRange); % Adjust color limits dynamically

    % Set axis limits and labels
    axis equal
    xlim([minX, maxX]);
    ylim([minY, maxY]);
    xlabel('Easting (km)'); 
    ylabel('Northing (km)');  
    set(gca, 'FontSize', 14, 'FontUnits', 'points');
    set(gca, 'XMinorTick', 'on');
    set(gca, 'YMinorTick', 'on');

    % Add the title
    title(elementTitles{i});

    nexttile;

    % Apply a mask to data to eliminate extrapolation
    Z = ApplyMask(Data(:,[2 1]), data(:,[1 2]), mean(datasimG(:, :, a),2), nx ,ny, 1);

    % Scatter plot with adjusted color mapping
    scatter(data(:, 2), data(:, 1), 20, Z, 'filled',"square");
    hold on
    plot(CoorHorne(2), CoorHorne(1), 'pk', 'MarkerSize', 10, 'MarkerFaceColor', 'k')
    hold on
    scatter(Data(:, 2), Data(:, 1), 20, 'k', 'filled', 'c');

    % Configure colorbar and colormap
    colorbar();
    colormap('jet'); % Change to 'balance' or 'redblue' if you prefer a zero-centered colormap
    clim(climRange); % Adjust color limits dynamically

    % Set axis limits and labels
    axis equal
    xlim([minX, maxX]);
    ylim([minY, maxY]);
    xlabel('Easting (km)');
    ylabel('Northing (km)');
    set(gca, 'FontSize', 14, 'FontUnits', 'points');
    set(gca, 'XMinorTick', 'on');
    set(gca, 'YMinorTick', 'on');

    % Add the title
    title(elementTitles{i});
end
%% Figure of std
a = 0;
figure(99)
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
elementTitles = { "log(As) (ppm)",  "log(Cd) (ppm)",};
for i = [1 2]
    a= a+1;
    nexttile;

    % Apply a mask to data to eliminate extrapolation
    Z = ApplyMask(Data(:,[2 1]), data(:,[1 2]), std(datasimNG(:, :, a),0,2), nx ,ny, 1);

    % Determine color limits dynamically to account for negative values
    if a ==1
        climRange = [0 1];
    else
        climRange = [0 2];
    end   

    % Scatter plot with adjusted color mapping
    scatter(data(:, 2), data(:, 1), 20, Z, 'filled',"square");
    hold on
    plot(CoorHorne(2), CoorHorne(1), 'pk', 'MarkerSize', 10, 'MarkerFaceColor', 'k')
    hold on
    scatter(Data(:, 2), Data(:, 1), 20, 'k', 'filled', 'c');

    % Configure colorbar and colormap
    colorbar();
    colormap('jet'); % Change to 'balance' or 'redblue' if you prefer a zero-centered colormap
    clim(climRange); % Adjust color limits dynamically

    % Set axis limits and labels
    axis equal
    xlim([minX, maxX]);
    ylim([minY, maxY]);
    xlabel('Easting (km)');
    ylabel('Northing (km)');
    set(gca, 'FontSize', 14, 'FontUnits', 'points');
    set(gca, 'XMinorTick', 'on');
    set(gca, 'YMinorTick', 'on');

    % Add the title
    title(elementTitles{i});

    nexttile;

    % Apply a mask to data to eliminate extrapolation
    Z = ApplyMask(Data(:,[2 1]), data(:,[1 2]), std(datasimG(:, :, a),0,2), nx ,ny, 1);

    % Scatter plot with adjusted color mapping
    scatter(data(:, 2), data(:, 1), 20, Z, 'filled',"square");
    hold on
    plot(CoorHorne(2), CoorHorne(1), 'pk', 'MarkerSize', 10, 'MarkerFaceColor', 'k')
    hold on
    scatter(Data(:, 2), Data(:, 1), 20, 'k', 'filled', 'c');

    % Configure colorbar and colormap
    colorbar();
    colormap('jet'); % Change to 'balance' or 'redblue' if you prefer a zero-centered colormap
    clim(climRange); % Adjust color limits dynamically

    % Set axis limits and labels
    axis equal
    xlim([minX, maxX]);
    ylim([minY, maxY]);
    xlabel('Easting (km)');
    ylabel('Northing (km)');
    set(gca, 'FontSize', 14, 'FontUnits', 'points');
    set(gca, 'XMinorTick', 'on');
    set(gca, 'YMinorTick', 'on');

    % Add the title
    title(elementTitles{i});
end

%% Figure of proba
a = 0;
figure(104)
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
elementTitles = { "P(As > 30 ppm)",  "P(Cd > 5 ppm)",};
for i = [1 2]
    a= a+1;
    nexttile;

    % Apply a mask to data to eliminate extrapolation
    if a ==1
        Z = ApplyMask(Data(:,[2 1]), data(:,[1 2]), mean(datasimNG_Org(:, :, a)>30,2), nx ,ny, 1);
    else
        Z = ApplyMask(Data(:,[2 1]), data(:,[1 2]), mean(datasimNG_Org(:, :, a)>5,2), nx ,ny, 1);
    end

    % Determine color limits dynamically to account for negative values
    climRange = [0 1];


    % Scatter plot with adjusted color mapping
    scatter(data(:, 2), data(:, 1), 20, Z, 'filled',"square");
    hold on
    plot(CoorHorne(2), CoorHorne(1), 'pk', 'MarkerSize', 10, 'MarkerFaceColor', 'k')
    hold on
    scatter(Data(:, 2), Data(:, 1), 20, 'k', 'filled', 'c');

    % Configure colorbar and colormap
    colorbar();
    colormap('jet'); % Change to 'balance' or 'redblue' if you prefer a zero-centered colormap
    clim(climRange); % Adjust color limits dynamically

    % Set axis limits and labels
    axis equal
    xlim([minX, maxX]);
    ylim([minY, maxY]);
    xlabel('Easting (km)');
    ylabel('Northing (km)');
    set(gca, 'FontSize', 14, 'FontUnits', 'points');
    set(gca, 'XMinorTick', 'on');
    set(gca, 'YMinorTick', 'on');

    % Add the title
    title(elementTitles{i});

    nexttile;

    % Apply a mask to data to eliminate extrapolation
    if a ==1
        Z = ApplyMask(Data(:,[2 1]), data(:,[1 2]), mean(datasimG_Org(:, :, a)>30,2), nx ,ny, 1);
    else
        Z = ApplyMask(Data(:,[2 1]), data(:,[1 2]), mean(datasimG_Org(:, :, a)>5,2), nx ,ny, 1);
    end

    % Determine color limits dynamically to account for negative values
    climRange = [0 1];

    % Scatter plot with adjusted color mapping
    scatter(data(:, 2), data(:, 1), 20, Z, 'filled',"square");
    hold on
    plot(CoorHorne(2), CoorHorne(1), 'pk', 'MarkerSize', 10, 'MarkerFaceColor', 'k')
    hold on
    scatter(Data(:, 2), Data(:, 1), 20, 'k', 'filled', 'c');

    % Configure colorbar and colormap
    colorbar();
    colormap('jet'); % Change to 'balance' or 'redblue' if you prefer a zero-centered colormap
    clim(climRange); % Adjust color limits dynamically

    % Set axis limits and labels
    axis equal
    xlim([minX, maxX]);
    ylim([minY, maxY]);
    xlabel('Easting (km)');
    ylabel('Northing (km)');
    set(gca, 'FontSize', 14, 'FontUnits', 'points');
    set(gca, 'XMinorTick', 'on');
    set(gca, 'YMinorTick', 'on');

    % Add the title
    title(elementTitles{i});
end

%% Cross-Validation
% Initialize error metrics
rmseNG = zeros(length([1 2]), 2); rmseG = zeros(length([1 2]), 2);
maeNG = zeros(length([1 2]), 2); maeG = zeros(length([1 2]), 2);
nrmseNG = zeros(length([1 2]), 2); nrmseG = zeros(length([1 2]), 2);

% Number of cross-validation steps
nstep = 10;

% Identify valid (non-NaN) samples
idx = all(~isnan(data(:, 2 + [1 2])), 2);
nn = sum(idx); % Number of valid samples

s = RandStream('mt19937ar','Seed',123);
cv = cvpartition(nn, 'KFold', nstep);  
RandStream.setGlobalStream(s);

nbsimul = 100;
for fold = 1:nstep

    a = 0;
    % Loop over each variable in [1 2]
    for i = [1 2]
        a = a + 1;
        idx = all(~isnan(data(:, 2 + [1 2])), 2);
        % Prepare training and validation datasets
        fullData = [x0(idx,1:2) data(idx, 2 + i)];

        trainIdx = training(cv, fold);
        testIdx = test(cv, fold);

        param.HD{a} = fullData(trainIdx, :); % Training set
        paramG.HD{a} = fullData(trainIdx, :); % Training set
        VC{a} = fullData(testIdx, :); % Validation set
        fullData(:,end) = exp(data(idx,end)*datastd(i) + datamean(i));
        VC_OD{a} = fullData(testIdx, :);
    end

    % Generate multivariate non-Gaussian random fields
    tic;
    nbfam = 21;
    [datasimNG_CV, Err] = GFFTMA_NG(model, c, nu, param, seed+fold*10, nbsimul, nbfam, nx, 1, ny, 1);
    tNG = toc;

    % Generate multivariate Gaussian random fields
    tic;
    datasimG_CV = GFFTMA_NG(model, c, nu, paramG, seed+fold*10, nbsimul, 1, nx, 1, ny, 1);
    tG = toc;

    % Apply inverse anamorphosis
    datasimG_CV_OD = datasimG_CV;
    datasimNG_CV_OD = datasimNG_CV;
    a = 0;
    for i=[1 2]
        a= a+1;
        datasimG_CV_OD(:, :, a) = exp(datasimG_CV(:,:,a)*datastd(i) + datamean(i));
        datasimNG_CV_OD(:, :, a) = exp(datasimNG_CV(:,:,a)*datastd(i) + datamean(i));
    end

    % Cross-Validation
    for a = 1:length([1 2])
        idx = sub2ind([nx, ny], VC{a}(:, 1), VC{a}(:, 2));

        % -------- Cas brut --------
        y_true = VC{a}(:, end);
        y_pred_G = datasimG_CV(idx, :, a);
        y_pred_NG = datasimNG_CV(idx, :, a);

        stats_G = regression_metrics(repmat(y_true, 1, nbsimul), y_pred_G);
        stats_NG = regression_metrics(repmat(y_true, 1, nbsimul), y_pred_NG);

        rmseG(a,1)    = rmseG(a,1)    + stats_G.RMSE / nstep;
        rmseNG(a,1)   = rmseNG(a,1)   + stats_NG.RMSE / nstep;
        maeG(a,1)     = maeG(a,1)     + stats_G.MAE / nstep;
        maeNG(a,1)    = maeNG(a,1)    + stats_NG.MAE / nstep;
        nrmseG(a,1)   = nrmseG(a,1)   + stats_G.nRMSE_std / nstep;
        nrmseNG(a,1)  = nrmseNG(a,1)  + stats_NG.nRMSE_std / nstep;

        % -------- Cas OD (transfo inverse) --------
        y_true_OD = VC_OD{a}(:, end);
        y_predG_OD = datasimG_CV_OD(idx, :, a);
        y_predNG_OD = datasimNG_CV_OD(idx, :, a);

        stats_G_OD = regression_metrics(repmat(y_true_OD, 1, nbsimul), y_predG_OD);
        stats_NG_OD = regression_metrics(repmat(y_true_OD, 1, nbsimul), y_predNG_OD);

        rmseG(a,2)    = rmseG(a,2)    + stats_G_OD.RMSE / nstep;
        rmseNG(a,2)   = rmseNG(a,2)   + stats_NG_OD.RMSE / nstep;
        maeG(a,2)     = maeG(a,2)     + stats_G_OD.MAE / nstep;
        maeNG(a,2)    = maeNG(a,2)    + stats_NG_OD.MAE / nstep;
        nrmseG(a,2)   = nrmseG(a,2)   + stats_G_OD.nRMSE_std / nstep;
        nrmseNG(a,2)  = nrmseNG(a,2)  + stats_NG_OD.nRMSE_std / nstep;
    end

end

%% Function %%

function [prediction_grid_xy, cols, rows] = make_grid(xmin, xmax, ymin, ymax, res)
    cols = round((xmax - xmin)/res);
    rows = round((ymax - ymin)/res);

    x = xmin + (1:cols)*res;
    y = ymin + (1:rows)*res;

    [xx, yy] = meshgrid(x, y);
    x = reshape(xx, [], 1);
    y = reshape(yy, [], 1);
    
    prediction_grid_xy = [x, y];
end

function [df_grid, grid_matrix, rows, cols] = grid_data(df, xx, yy, zz, res)
    % Rename columns
    X = df.(xx);
    Y = df.(yy);
    Z = df.(zz);

    % Bounding box
    xmin = min(X)-res/2; xmax = max(X)+res/2;
    ymin = min(Y)-res/2; ymax = max(Y)+res/2;

    % Generate grid
    [grid_coord, cols, rows] = make_grid(xmin, xmax, ymin, ymax, res);

    % Convert to arrays
    np_data = [X, Y, Z];
    np_resize = np_data;
    origin = [xmin, ymin];
    resolution = [res, res];

    % Rescale positions to grid indices
    np_resize(:, 1:2) = round((np_resize(:, 1:2) - origin) ./ resolution);

    % Init
    grid_sum = zeros(rows, cols);
    grid_count = zeros(rows, cols);

    % Accumulate
    for i = 1:size(np_data, 1)
        xindex = int32(np_resize(i, 2)) + 1;
        yindex = int32(np_resize(i, 1)) + 1;

        if (xindex > rows || yindex > cols || xindex < 1 || yindex < 1)
            continue;
        end

        grid_sum(xindex, yindex) = grid_sum(xindex, yindex) + np_data(i, 3);
        grid_count(xindex, yindex) = grid_count(xindex, yindex) + 1;
    end

    % Compute average
    grid_matrix = grid_sum ./ grid_count;
    grid_matrix(grid_count == 0) = NaN;

    % Flatten everything
    grid_array     = reshape(grid_matrix, [], 1);
    grid_sum_flat  = reshape(grid_sum, [], 1);
    grid_count_flat = reshape(grid_count, [], 1);

    % Make sure grid_coord is reshaped to match
    if size(grid_coord,1) ~= length(grid_array)
        grid_coord = reshape(grid_coord, [], 2);  % safety
    end

    % Create table
    df_grid = table(grid_coord(:,1), grid_coord(:,2), ...
                    grid_sum_flat, grid_count_flat, grid_array, ...
                    'VariableNames', {'X', 'Y', 'Sum', 'Count', 'Z'});

    % Flip matrix vertically to match display
    grid_matrix = grid_matrix;
end

function metrics = regression_metrics(y_true, y_pred)
    % Vérification de la taille des données
    [n, m] = size(y_true); % n = nombre de points de données, m = nombre de simulations
    if size(y_pred) ~= size(y_true)
        error('y_true et y_pred doivent avoir la même taille.');
    end
    
    % Calcul des résidus
    residuals = y_pred - y_true;
    
    % Erreur absolue moyenne (MAE) sur chaque simulation
    mae = mean(abs(residuals), 1); % Moyenne des erreurs absolues par colonne (simulation)
    
    % Erreur quadratique moyenne (RMSE) sur chaque simulation
    rmse = sqrt(mean(residuals.^2, 1)); % Moyenne des erreurs quadratiques par colonne
    
    % nRMSE normalisé par la moyenne
    nrmse_mean = rmse / mean(y_true, 1); % Normalisé par la moyenne (par colonne)
    
    % nRMSE normalisé par l'étendue (max - min)
    nrmse_range = rmse / (max(y_true, [], 1) - min(y_true, [], 1)); % Normalisé par l'étendue
    
    % nRMSE normalisé par l'écart-type
    nrmse_std = sqrt(mean((residuals ./ std(y_pred, 0, 1)).^2, 1)); % Normalisé par l'écart-type
    
    % MAPE calculation (Mean Absolute Percentage Error)
    numerator = abs(residuals);
    denominator =  abs(y_true);

    % Éviter la division par zéro
    zero_mask = denominator == 0;
    numerator(zero_mask) = 0;
    denominator(zero_mask) = 1;

    mape = 100 * mean(numerator ./ denominator, 1); % Moyenne par colonne


    % SMAPE (Symmetric Mean Absolute Percentage Error)
    numerator = abs(residuals);
    denominator = abs(y_pred) + abs(y_true);
    
    % Éviter la division par zéro
    zero_mask = denominator == 0;
    numerator(zero_mask) = 0;
    denominator(zero_mask) = 1;
    
    smape = 100 * mean(numerator ./ denominator, 1); % Moyenne par colonne

    % R² (Coefficient de détermination) sur chaque simulation
    ss_res = sum(residuals.^2, 1); % Somme des carrés des résidus
    ss_tot = sum((y_true - mean(y_true, 1)).^2, 1); % Somme des carrés totaux
    r2 = 1 - (ss_res ./ ss_tot); % R² par simulation

    % AIC / BIC
    k = 1; % Nombre de paramètres estimés (ajustement simple ici)
    logL = -n/2 * log(2*pi) - n/2 * log(sum(residuals.^2, 1) / n) - n/2;
    AIC = -2*logL + 2*k; % AIC
    BIC = -2*logL + log(n)*k; % BIC

    % Résultat final : Moyenne de chaque métrique sur les simulations
    metrics = struct(...
        'MAE', mean(mae), ...
        'RMSE', mean(rmse), ...
        'nRMSE_mean', mean(nrmse_mean), ...
        'nRMSE_range', mean(nrmse_range), ...
        'nRMSE_std', mean(nrmse_std), ...
        'MAPE', mean(mape), ...
        'SMAPE', mean(smape), ...
        'R2', mean(r2), ...
        'AIC', mean(AIC), ...
        'BIC', mean(BIC) ...
    );
end



