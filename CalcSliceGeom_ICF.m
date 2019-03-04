function [ predP, predR ] = CalcSliceGeom_ICF( I, z)
%CALCPVE_VP Summary of this function goes here
%   Detailed explanation goes here


% Scale factors for tilted cylinder
zXY = atan(z(2)./z(1));
if isnan(zXY)
    zXY = 0;
end
xScale = sqrt((sin(zXY)./z(3)).^2 + (cos(zXY)).^2);
yScale = sqrt((cos(zXY)./z(3)).^2 + (sin(zXY)).^2);

% x dimension
mass = cumsum(sum(I,1));
massFrac = mass./sum(I(:));
midPoint = sum(I(:))/2;
xBisect1 = find(mass<midPoint,1,'last');
xBisect2 = find(mass>midPoint,1,'first');

t1x = TSinT_ICF(massFrac(xBisect1));
t2x = TSinT_ICF(1-massFrac(xBisect2));
rx = (xBisect2-xBisect1)./(cos(t1x/2)+cos(t2x/2));

% y dimension
mass = cumsum(sum(I,2));
massFrac = mass./sum(I(:));
midPoint = sum(I(:))/2;
yBisect1 = find(mass<midPoint,1,'last');
yBisect2 = find(mass>midPoint,1,'first');

t1y = TSinT_ICF(massFrac(yBisect1));
t2y = TSinT_ICF(1-massFrac(yBisect2));
ry = (yBisect2-yBisect1)./(cos(t1y/2)+cos(t2y/2));

x1 = xBisect1 + rx*cos(t1x/2);
x2 = xBisect2 - rx*cos(t2x/2);
x = 0.5 * (x1+ x2 + 2);

y1 = yBisect1 + ry*cos(t1y/2);
y2 = yBisect2 - ry*cos(t2y/2);
y = 0.5 * (y1 + y2 + 2);

%r = 0.5*(rx.*xScale+ry.*yScale);

predP = [y x]-1; %correct for 0,0 origin not 1,1
%predR = 0.5*(ry./(1+cos(zXY).^2) + rx./(1+sin(zXY).^2));%.*(sin(zXY).^2);%./(1+(z(3).^2).*(sin(zXY.*2).^2));
predR = 0.5*(ry./xScale+rx./yScale);
end
