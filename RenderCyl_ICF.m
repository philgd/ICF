function [ partialVol ] = RenderCyl_ICF( dims, cylO, cylR, cylZ, quick )
%RENDER_VP Summary of this function goes here
%   Detailed explanation goes here

rng('shuffle')

nPoints = 400;

[x,y,z] = ind2sub(dims,1:prod(dims));
x = repmat(x-1,[nPoints 1])+rand([nPoints numel(x)]);
y = repmat(y-1,[nPoints 1])+rand([nPoints numel(y)]);
z = repmat(z-1,[nPoints 1])+rand([nPoints numel(z)]);
[ pOut ] = CartToCylindrical_ICF( [x(:) y(:) z(:)], cylO, cylZ);
rad = reshape(pOut(:,1)<cylR,[nPoints prod(dims)]);
partialVol = reshape(sum(rad,1)./nPoints,dims);

% if nargin<5
%     quick = true;
% end
% 
% if numel(dims) == 2
%     dims = [dims 1];
% end
% 
% z = z(:);
% p = p(:);
% 
% kStart = 1;
% kEnd = dims(3);
% 
% % pVol == 1 if all corners are inside circle, no need to test further
% % Within 1.5 pixels (+diags) are partial, otherwise zero
% if quick
%     partialVol = zeros(dims);
%     for k=kStart:kEnd
%         pi = p(1:2) + (k-p(3)-1)*z([1 2])/z(3).*[1 -1]';
%         [ estMap ] = EstIntCyl_PV( dims([1 2]), pi, r, z);
%         partialVol(:,:,k) = estMap;
%     end
%     
%     pMask = partialVol == 1;
%     partialVol = 0.5*imdilate(partialVol,ball(1.5));
%     partialVol(pMask) = 1;
% else
%     partialVol = 0.5*ones(dims);
% end
% 
% % Still need to test partials
% for k=kStart:kEnd
%     pi = p(1:2) + (k-p(3)-1)*z([1 2])/z(3).*[1 -1]';
%     
%     pSlice = squeeze(partialVol(:,:,k));
%     pMask = squeeze(partialVol(:,:,k)>0 & partialVol(:,:,k)<1);
%     pSlice(pMask) = IntCyl_PV(pMask, pi, r, z);
%     
%     partialVol(:,:,k) = pSlice;
% end

end

