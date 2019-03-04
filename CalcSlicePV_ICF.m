function [ partialVol ] = CalcSlicePV_ICF( dims, cylO, cylR, cylZ)
%RENDER_VP Summary of this function goes here
%   Detailed explanation goes here

if size(cylZ,1)==1
    cylZ = repmat(cylZ,[size(cylO,1) 1]);
end


rng('shuffle')

nPoints = 500;

[x,y,z] = ind2sub(dims,1:prod(dims));
sliceInd = repmat(z,[nPoints 1]);
x = repmat(x-1,[nPoints 1])+rand([nPoints numel(x)]);
y = repmat(y-1,[nPoints 1])+rand([nPoints numel(y)]);
z = repmat(z-1,[nPoints 1])+rand([nPoints numel(z)]);

partialVol = zeros(dims);
for i=1:dims(3)
    mask = sliceInd==i;
    [ pOut ] = CartToCylindrical_ICF( [x(mask) y(mask) z(mask)], cylO(i,:), cylZ(i,:));
    rad = reshape(pOut(:,1)<cylR(i),[nPoints prod(dims([1 2]))]);
    partialVol(:,:,i) = reshape(sum(rad,1)./nPoints,dims([1 2]));

end
