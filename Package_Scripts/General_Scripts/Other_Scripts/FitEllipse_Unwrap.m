function el = FitEllipse_Unwrap(XY,elCenter)
%%
% Name: DirectEllipseFit_Unwrap
% Author: Daniel O'Hara
% Data: 03/03/2021 (mm/dd/yyyy)
% Description: Wrapper script for EllipseDirectFit code to calculate the
% best-fitting ellipse to the data points given in XY.
%
% Input:
%   XY: 2xN matrix of point x- and y-coordinates (first and second columns,
%       respectively).
%   elCenter: Ellipse center coordinates. Empty unless it is being used to
%       determine the ellipse of non-closed contours 
%
% Output: 
%   el: Class of ellipse parameters, including:
%       x0: Ellipse centroid x-value.
%       y0: Ellipse centroid y-value.
%       longAxis: Major axis radius.
%       shortAxis: Minor axis radius.
%       phi: Major axis orientation in azimuth.
%       mathPhi: Major axis orientation in mathematical radians.

try
    if nargin == 1 || isempty(elCenter)
        a = fitellipse(XY(:,1),XY(:,2));

        el.x0 = a(1);
        el.y0 = a(2);
        el.longAxis = a(3);
        el.shortAxis = a(4);
        el.mathPhi = rad2deg(a(5));
    else
        a = fit_ellipse_zero_R1(XY(:,1)-elCenter(1),XY(:,2)-elCenter(2));

        el.x0 = elCenter(1);
        el.y0 = elCenter(2);
        el.longAxis = a.long_radius;
        el.shortAxis = a.short_radius;
        el.mathPhi = 180-a.angleToX;
    end
    
    el.phi = 450 - el.mathPhi;

    while el.phi >180
        el.phi = el.phi - 180;
    end

    if el.longAxis < el.shortAxis
        tt = el.longAxis;
        el.longAxis = el.shortAxis;
        el.shortAxis = tt;
    end
    
catch er
    warning('No Ellipse Found')

    el.x0 = NaN;
    el.y0 = NaN;
    el.longAxis = NaN;
    el.shortAxis = NaN;
    el.mathPhi = NaN;
    el.phi = NaN;
end