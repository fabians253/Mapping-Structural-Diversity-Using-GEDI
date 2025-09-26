
%% Mapping-Structural-Diversity-Using-GEDI
% 
% main author: Fabian D. Schneider
% 
% This is the main Matlab script for mapping structural diversity using GEDI.
% 
% Please read and reference (cite) the following scientific paper when using this code:
% 
% Fabian D. Schneider*, Morgan Dean, Elsa M. Ordway, Moses B. Libalah, & Antonio A. Ferraz. Mapping the structural diversity of Central African and Western US forests using GEDI. In Review at Remote Sensing of Environment.
% *fabian.schneider@bio.au.dk; Section for Ecoinformatics & Biodiversity, Department of Biology, Aarhus University, Ny Munkegade 114, DK-8000 Aarhus, Denmark

%% Load and prepare Datasets

% load the GEDI data table with quality controlled GEDI L2A and L2B data
% downloaded from Google Earth Engine
load('data/KingsCanyon_GEDIL2AB_20190325_20230301.mat');

% read GEDI reference grid, for example from gridded biomass product
[gedi_agb, R_gedi_agb] = readgeoraster( 'data/KingsCanyon_GEDI04_B_MW019MW223_02_002_02_R01000M_MU.tif' );

% load GEDI trait ranges to be used for diversity mapping
load( 'data/GEDI_traits_ranges.mat' );

% define a 3D sample grid to calculate probability density in the 3D
% feature space, could be changed to any dimenion of interest

% 3D: 8000 samples
vec1 = 0:0.05:1;
xi3d = combvec( vec1, vec1, vec1 )';

% save the GEDI traits to be included in the diversity calculation
data_tab(:,1) = gedi_points_kingsCanyon.rh98;
data_tab(:,2) = gedi_points_kingsCanyon.cover;
data_tab(:,3) = gedi_points_kingsCanyon.fhd;

% save the GEDI coordiantes
x = gedi_points_kingsCanyon.x;
y = gedi_points_kingsCanyon.y;

% normalize traits and remove outliers if needed
minVals = [gedi_traits_ranges.rh98(1) gedi_traits_ranges.cover(1) gedi_traits_ranges.fhd(1)];
maxVals = [gedi_traits_ranges.rh98(2) gedi_traits_ranges.cover(2) gedi_traits_ranges.fhd(2)];
data_tab = scaleZeroOne_cols_absolute( data_tab, minVals, maxVals ); % normalize based on predefined trait ranges
% data_tab = scaleNoOutliers_cols( data_tab, 0.05 ); % remove a percentile of outliers such as 0.05%

%% Map Diversity

% define a probability density threshold for deriving richness
probThr = 0.2;

% define minimum number of points to calculate diversity
minPoints = 10;

% calculate diversity metrics for a particular resolution res in
% kilometers: here it is calculated for 1 and 5 km; since the test data is
% 5x5 km, the 5 km scale will end up in 1 pixel;
% Diversity metrics are based on the concept of functional richness (fric),
% evenness (feve) and divergenxe (fdiv)
for res = [1 5]

temp = R_gedi_agb.CellExtentInWorldX;
pixelSize = temp*res;

% get GEDI grid and georeference:
gedi_x_vec = R_gedi_agb.XWorldLimits(1):pixelSize:R_gedi_agb.XWorldLimits(2);
gedi_y_vec = fliplr( R_gedi_agb.YWorldLimits(1):pixelSize:R_gedi_agb.YWorldLimits(2) );

nrRows = length( gedi_y_vec ) - 1;
nrCols = length( gedi_x_vec ) - 1;

% initialize outputs
gedi_fric = zeros( nrRows, nrCols ) * NaN;
gedi_feve = zeros( nrRows, nrCols ) * NaN;
gedi_fdiv = zeros( nrRows, nrCols ) * NaN;
nrShots = zeros( nrRows, nrCols ) * NaN;

% this code is parallelized with parfor
t1 = tic;
for row = 1:nrRows
    disp( ['Row: ' num2str(row) ' / ' num2str( nrRows )] );
    t2 = tic;
    indRows = y < gedi_y_vec( row ) & y > gedi_y_vec( row+1 );
    %temp = data_tab( indRows, : );
    temp = parallel.pool.Constant( data_tab( indRows, : ) );
    x_temp = x( indRows );
    parfor col = 1:nrCols
        indCols = x_temp > gedi_x_vec( col ) & x_temp < gedi_x_vec( col+1 );
        X = temp.Value( indCols, : );

        % calculate functional richness, evenness and divergence
        [gedi_fric( row, col ), gedi_feve( row, col ), gedi_fdiv( row, col )] = getFRicFEveFDiv_PDFadapt( X, xi3d, probThr, minPoints );

        % save number of pixels used
        nans = any( isnan(X), 2 );
        X = unique( X( ~nans, : ), 'rows' );
        nrShots( row, col ) = size(X,1);
    end
    timeByRowSec(row) = toc(t2);
end
totalTimeMin = toc(t1)/60;

% show figure
figure;

subplot(2,2,1);
imagesc( gedi_fric ); axis image; colorbar
title('Richness');

subplot(2,2,2);
imagesc( gedi_feve ); axis image; colorbar
title('Evenness');

subplot(2,2,3);
imagesc( gedi_fdiv ); axis image; colorbar
title('Divergence');

subplot(2,2,4);
imagesc( nrShots ); axis image; colorbar
title('Nr Shots');

% save output
save( ['GEDI_Diversity_' num2str(res) 'km_SchneiderEtAl.mat'], 'gedi_fric', 'gedi_fdiv', 'gedi_feve', 'nrShots', 'pixelSize', "gedi_x_vec", "gedi_y_vec" );

end

%% Export to GeoTIFF

for spScale = [1 5]

load( ['GEDI_Diversity_' num2str(spScale) 'km_SchneiderEtAl.mat'] );

data = cat( 3, gedi_fric, gedi_feve, gedi_fdiv );
xEdges = gedi_x_vec([1 end]);
yEdges = gedi_y_vec([end 1]);
cellSize = pixelSize;
coordRefSysCode = 'EPSG:6933';
filePathOut = 'geotiff_out/';
fileName = ['GEDI_FRic_FEve_FDiv_' num2str(spScale) 'km_SchneiderEtAl.tif'];
exportGeotiff( data, xEdges, yEdges, cellSize, coordRefSysCode, filePathOut, fileName );

end
