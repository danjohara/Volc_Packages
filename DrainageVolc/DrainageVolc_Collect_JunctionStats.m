function J_X_Y_C_Z_D_A = DrainageVolc_Collect_JunctionStats(DEM,DVD,As,CJ,X,Y)
% Name: collectJunctionStats
% Author: Daniel O'Hara
% Date: 02/16/2021 (mm/dd/yyyy)
% Description: Script to determine juction statistics for drainage divides.
%
% Input:
%   DEM: GRIDobj of elevations.
%   DVD: DIVIDEobj of draiage divides.
%   A: List of asymmetry indices.
%   CJ: List of divide junction connectivity values.
%   X: X-coordinates of divide junctions.
%   Y: Y-coordinates of divide junctions.
%
% Output:
%   J_X_Y_C_Z_D_A: Matrix divide junction values. Columns correspond to
%       junction x-coordinate, y-coordinate, junction connectivity,
%       elevation, divide distance, and asymmetry index.

J_X_Y_C_Z_D_A = zeros([length(CJ),6]);
J_X_Y_C_Z_D_A(:,1:3) = [X,Y,CJ];

[zz,xx,yy] = GRIDobj2mat(DEM);
[xx,yy] = meshgrid(xx,yy);
zzs = interp2(xx,yy,zz,X,Y);

for i = 1:length(CJ)
    ind = coord2ind(DVD,X(i),Y(i));
    z = zzs(i);
    listI = find(DVD.IX==ind,1);
    dist = DVD.distance(listI);
    asy = As(listI);
    
    J_X_Y_C_Z_D_A(i,4) = z;
    J_X_Y_C_Z_D_A(i,5) = dist;
    J_X_Y_C_Z_D_A(i,6) = asy;
end