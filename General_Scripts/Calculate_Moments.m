function [ x,y,z,Gm,PDFm,centerOfMass,skewness,kurtosis ] = Calculate_Moments(x,y,z)
% Name: calcMoments3
% Author: Daniel O'Hara
% Date: 02/04/2021 (mm/dd/yyyy) - original script created 2019.
% Description: Script to calculate the mathematical moments of shape.
%
% Input:
%   x: Vector of x coordinates.
%   y: Vector of y coordinates.
%   z: Grid of Z values.
%
% Output:
%   x: Expanded vector of x coordinates.
%   y: Expanded vector of y coordinates.
%   z: Expanded grid of Z values.
%   PDFm: Matrix of moments, where each row/column corresponds to the
%       0-based order of the moment (e.g., first row/column is the 0th 
%       order). X-based moments are columns, Y-based moments are rows.
%   centerOfMass: X-Y coordinates of the volumetric center of mass
%       (centroid).
%   skewness: Values of mathematical skewness in X and Y directions.
%       Negative values imply skewness in the West/South direction,
%       positive values imply skewness in the East/North direction.
%   kurtosis: Values of mathematical kurtosis (i.e. 'tailedness' of the
%       shape in the X and Y directions.


dx = abs(diff(x(1,1:2)));
dy = abs(diff(y(1,1:2)));
oldX = x;
oldY = y;

x = [min(x)-dx,x,max(x)+dx];
y = [min(y)-dx,y,max(y)+dy];

z0s = zeros(size(oldX));
z0s2 = [0;zeros(size(oldY))';0];

z = [z0s2,[z0s;z;z0s],z0s2];

sx = size(x);
sy = size(y);

xG = repmat(x,[sy(1,2) 1]);
yG = repmat(y',[1,sx(1,2)]);

%Geometric moments
Gm = zeros(3,3);

for i = 1:3
    for j = 1:3
        zx = trapz(yG.^(j-1).*xG.^(i-1).*z,2);
        zy = trapz(zx,1);
        Gm(j,i) = zy;
    end
end

centerOfMass = [Gm(1,2)/Gm(1,1), Gm(2,1)/Gm(1,1)];

%PDF moments
PDFm = zeros(5,5);

% Centralize to the CoM
xG = xG - centerOfMass(1);
yG = yG - centerOfMass(2);

% z = z./Gm(1,1);

for i = 1:5
    for j = 1:5
        zx = trapz(yG.^(j-1).*xG.^(i-1).*z,2);
        zy = trapz(zx,1);
        PDFm(j,i) = zy;
    end
end

% lll = ones(5,5);
% lll(1,2:5) = sqrt(PDFm(1,3));
% lll(2:5,1) = sqrt(PDFm(3,1));
% mmm = zeros(5,5);
% mmm(1,2:5) = 1:4;
% mmm(2:5,1) = 1:4;
% lmlm = lll.^mmm;
% 
% ppp = PDFm./lmlm;
% 
% skewness = [ppp(1,4),ppp(4,1)];
% kurtosis = [ppp(1,5),ppp(5,1)];

% Normalize to the 0th moment.
PDFm = PDFm./PDFm(1,1);

skewness = [PDFm(1,4)/sqrt(PDFm(1,3)^3),PDFm(4,1)/sqrt(PDFm(3,1)^3)];
kurtosis = [PDFm(1,5)/PDFm(1,3)^2,PDFm(5,1)/PDFm(3,1)^2];
end

