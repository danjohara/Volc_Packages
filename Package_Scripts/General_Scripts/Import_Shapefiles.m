function [boundaryXY,craterXY,maskXY] = Import_Shapefiles(boundaryXY,craterXY,maskXY,verbose)
%%
% Name: Import_Shapefiles
% Author: Daniel O'Hara
% Data: 05/22/2024 (mm/dd/yyyy)
% Description: Simple function to convert boundary, crater, and mask files
%   into XY data usable by the Volc packages.
%
% Input:
%   boundaryXY: Boundary information. If given as a shapefile, code will
%       import the data with shaperead; otherwise, the code will assume it
%       is already in its proper format as an Nx2 matrix (first row is X,
%       second row is Y). In all cases, code will try to convert data to
%       UTM format. If boundaryXY is given as a matrix already and in
%       geographic coodinates, code assumes longitude is still in the first
%       row.
%   craterXY: Crater information. File is treated the same as boundaryXY.
%       Code accepts multiple craters; either as multiple entries in the
%       shapefile, or as a cell array.
%   maskXY: Mask information. File is treated the same as boundaryXY.
%       Code accepts multiple mask regions; either as multiple entries in 
%       the shapefile, or as a cell array.
%   verbose: Flag for script progression output.
%
% Output: 
%   boundaryXY: Nx2 matrix of boundary XY values.
%   craterXY: Cell array of Nx2 matrices of crater XY values.
%   maskXY: Cell array of Nx2 matrices of mask XY values.

%% Import boundary file
    if ~isempty(boundaryXY)
        if verbose > 0
            disp('Importing boundary data...')
        end

        % If boundary is given as shapefile, convert to an array.
        if ischar(boundaryXY)
            Sh = shaperead(boundaryXY);
            tx = Sh.X;
            ty = Sh.Y;
            if isnan(tx(end))
                tx(end) = [];
                ty(end) = [];
            end
    
            boundaryXY = [tx',ty'];
        end 
    
        if nansum((abs(boundaryXY(:,1))>180)*1) == 0 || nansum((abs(boundaryXY(:,2))>90)*1) == 0
            try
                [tx,ty,zone] = ll2utm(boundaryXY(:,2),boundaryXY(:,1));
                if length(unique(zone)) > 1
                    [tx,ty,~] = ll2utm(boundaryXY(:,2),boundaryXY(:,1),mode(zone));
                end

                boundaryXY = [tx,ty];
            catch
                warning('Unable to convert boundary from Lat/Lon, assuming already in UTM')
            end
        end
    end
    
%% Import crater file
    if ~isempty(craterXY)
        if verbose > 0
            disp('Importing crater data...')
        end

        % If crater is given as shapefile, convert to cell array.
        if ischar(craterXY)
            Sh = shaperead(craterXY);
            craterXY = cell(length(Sh),1);
            for i = 1:size(Sh,1)
                tx = Sh(i).X;
                ty = Sh(i).Y;
                if isnan(tx(end))
                    tx(end) = [];
                    ty(end) = [];
                end
    
                craterXY{i} = [tx',ty'];
            end
        elseif ~isempty(craterXY) && ~iscell(craterXY)
            craterXY = {craterXY};
        end
    
        for i = 1:length(craterXY)
            if nansum((abs(craterXY{i}(:,1))>180)*1) == 0 || nansum((abs(craterXY{i}(:,2))>90)*1) == 0
                try
                    [tx,ty,zone] = ll2utm(craterXY{i}(:,2),craterXY{i}(:,1));
                    if length(unique(zone)) > 1
                        [tx,ty,~] = ll2utm(craterXY{i}(:,2),craterXY{i}(:,1),mode(zone));
                    end

                    craterXY{i} = [tx,ty];
                catch
                    warning('Unable to convert crater from Lat/Lon, assuming already in UTM')
                end
            end
        end
    end
    
%% Import mask file
    if ~isempty(maskXY)
        if verbose > 0
            disp('Importing mask data...')
        end

        % If mask is given as shapefile, convert to an array.
        if ischar(maskXY)
            Sh = shaperead(maskXY);
            maskXY = cell(length(Sh),1);
            for i = 1:size(Sh,1)
                tx = Sh(i).X;
                ty = Sh(i).Y;
                if isnan(tx(end))
                    tx(end) = [];
                    ty(end) = [];
                end
    
                maskXY{i} = [tx',ty'];
            end
        elseif ~isempty(maskXY) && ~iscell(maskXY)
            maskXY = {maskXY};
        end
    
        for i = 1:length(maskXY)
            if nansum((abs(maskXY{i}(:,1))>180)*1) == 0 || nansum((abs(maskXY{i}(:,2))>90)*1) == 0
                try
                    [tx,ty,zone] = ll2utm(maskXY{i}(:,2),maskXY{i}(:,1));
                    if length(unique(zone)) > 1
                        [tx,ty,~] = ll2utm(maskXY{i}(:,2),maskXY{i}(:,1),mode(zone));
                    end
                    maskXY{i} = [tx,ty];
                catch
                    warning('Unable to convert mask from Lat/Lon, assuming already in UTM')
                end
            end
        end
    end
end