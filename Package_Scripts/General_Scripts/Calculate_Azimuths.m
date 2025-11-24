function az = Calculate_Azimuths(oXY,XY)
%%
% Name: Calculate_Azimuths
% Initial Date: 3/04/2025 (mm/dd/yyyy)
% Author: Daniel O'Hara
% Description: Simple script to calculate azimuthal directions of points 
%   from an origin.
%
% Input: 
%   oXY: 1x2 array of X-Y coordinates of the origin.
%   XY: Mx2 matrix of X-Y coordinates of points to calculate azimuth from 
%       origin.
%   
% Output:
%   az: Mx1 array of azimuths.

dXs = XY(:,1)-oXY(1);
dYs = XY(:,2)-oXY(2);

az = atan2(dYs,dXs);
az = rad2deg(az);
az(az<0) = az(az<0) + 360;
az = 450 - az;
az(az>360) = az(az>360) - 360;
end