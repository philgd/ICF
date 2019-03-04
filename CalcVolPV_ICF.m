function [ partialVol ] = CalcVolPV_ICF( dims, cylO, cylR, cylZ)
%RENDER_VP Summary of this function goes here
%   Detailed explanation goes here

rng('shuffle')

nPoints = 500;

[x,y,z] = ind2sub(dims,1:prod(dims));
x = repmat(x-1,[nPoints 1])+rand([nPoints numel(x)]);
y = repmat(y-1,[nPoints 1])+rand([nPoints numel(y)]);
z = repmat(z-1,[nPoints 1])+rand([nPoints numel(z)]);
[ pOut ] = CartToCylindrical_ICF( [x(:) y(:) z(:)], cylO, cylZ);
rad = reshape(pOut(:,1)<cylR,[nPoints prod(dims)]);
partialVol = reshape(sum(rad,1)./nPoints,dims);

end

