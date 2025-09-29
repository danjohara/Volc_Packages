function [DEM0,X0,Y0,Z0,DEM,X,Y,Z] = Import_DEM(tifFile,dx,fillDEM,boundaryXY,craterXY,maskXY,verbose)
%%
% Name: Import_DEM
% Author: Daniel O'Hara
% Data: 05/22/2024 (mm/dd/yyyy)
% Description: Function to import raster data, fill sinks (if appropriate),
%   and cut data for use in the Volc packages.
%
% Input:
%   tifFile: Raster data. Can either be a GRIDobj variable or a folder path
%       to a .tif file. In all cases, script will attempt to reproject the
%       grid to UTM coordinates, and resample to given cell size.
%   dx: Grid resolution.
%   fillDEM: Flag for whether the raster should be filled (if being used
%       for drainage analysis, such as DrainageVolc).
%   boundaryXY: Nx2 matrix of boundary XY values.
%   craterXY: Cell array of Nx2 matrices of crater XY values.
%   maskXY: Cell array of Nx2 matrices of mask XY values.
%   verbose: Flag for script progression output.
%
% Output: 
%   DEM0: Resampled GRIDobj of initial raster data.
%   X, Y, Z: X-Y coordinates and elevations of intial raster data (with
%       sinks filled if appropriate).
%   DEM: Cut and masked GRIDobj of raster data.
%   cutX, cutY, cutZ: X-Y coordinates and elevations of cut/masked raster
%       data.

%% Import DEM
    if verbose > 0
        disp('Loading DEM...')
    end

    if ischar(tifFile)
        warning('off','all')
        DEM0 = GRIDobj(tifFile);
        warning('on','all')
    
        DEM0.Z(DEM0.Z<-20000) = NaN;
    else
        DEM0 = tifFile;
    end

    V = 2;
    try 
        DEM0.refmat;
    catch
        V = 3;
    end

    if V==2 && isempty(DEM0.georef)
        try
            DEM0 = reproject2utm(DEM0,dx);
        catch
            warning('Cannot reproject DEM to UTM, assuming map is already in coordinates.')
            DEM0 = resample(DEM0,dx);
        end
    elseif V==3 && isempty(DEM0.georef.ProjectedCRS)
        try
            DEM0 = reproject2utm(DEM0,dx);
        catch
            warning('Cannot reproject DEM to UTM, assuming map is already in coordinates.')
            DEM0 = resample(DEM0,dx);
        end
    else
        DEM0 = resample(DEM0,dx);
    end

    % TopoToolbox 3 resampling does not always create uniformly-spaced
    % grids, which 
    if V==3 && (DEM0.wf(1,1) ~= DEM0.wf(2,2))
        warning('Resample created non-uniform coordinates; re-interpolating.')
        tmpDEM = DEM0;
        [Z,x0,y0] = GRIDobj2mat(DEM0);
        x1 = min(x0):dx:max(x0);
        y1 = min(y0):dx:max(y0);

        [X0,Y0] = meshgrid(x0,y0);
        [X1,Y1] = meshgrid(x1,y1);
        scInt = scatteredInterpolant(double(X0(:)),double(Y0(:)),double(Z(:)),'natural','none');
        DEM0 = GRIDobj(X1,Y1,scInt(X1,Y1));
        DEM0.georef.ProjectedCRS = tmpDEM.georef.ProjectedCRS;
    end

%% Fill DEM
    DEM = DEM0;
    if fillDEM
        if verbose > 0
            disp('Filling DEM...')
        end
        if isempty(craterXY)
            DEM = fillsinks(DEM);
        else
            % Ignore crater regions when filling sinks.
            [~,tx,ty] = GRIDobj2mat(DEM);
            [tX,tY] = meshgrid(tx,ty);
            tS = DEM;
            tS.Z = zeros(size(tS.Z));
            for i = 1:length(craterXY)
                inP = inpolygon(tX,tY,craterXY{i}(:,1),craterXY{i}(:,2));
                tS.Z(inP==1) = 1;
            end
            DEM = fillsinks(DEM,tS);
        end
    end
    
%% Get big DEM values
    [Z0,X0,Y0] = GRIDobj2mat(DEM);
    [XG0,YG0] = meshgrid(X0,Y0);
    
%% Apply mask and cut around boundary
    if verbose > 0
        disp('Masking and Clipping Topography...')
    end
    newZ = Z0;
    
    if ~isempty(boundaryXY)
        inP = inpolygon(XG0,YG0,boundaryXY(:,1),boundaryXY(:,2));
        newZ(~inP) = NaN;
    end
    
    for i = 1:length(maskXY)
        inP = inpolygon(XG0,YG0,maskXY{i}(:,1),maskXY{i}(:,2));
        newZ(inP) = NaN;
    end
    
    DEM = GRIDobj(XG0,YG0,newZ);
    DEM.georef = DEM0.georef;

    if isfield(DEM0,'refmat')
        DEM.refmat = DEM0.refmat;
    elseif isfield(DEM0,'wf')
        DEM.wf = DEM0.wf;
    end
    
    DEM = crop(DEM);
    [Z,X,Y] = GRIDobj2mat(DEM);

%% Set everything as double
    X0 = double(X0);
    Y0 = double(Y0);
    Z0 = double(Z0);
    X = double(X);
    Y = double(Y);
    Z = double(Z);
end

