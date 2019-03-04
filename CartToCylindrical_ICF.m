function [ pOut ] = CartToCylindrical_ICF( p, cylO, cylZ)
%CARTTOCYLINDRICAL Summary of this function goes here
%   Detailed explanation goes here

pOut = zeros(size(p));
pTemp = zeros(size(p));

% Correct for cylinder origin
p = p-repmat(cylO,[size(p,1) 1]);

% Magnitude of vector parallel to the cylinder (height)
paraM = p*cylZ';
pOut(:,3) = paraM;

% Magnitude of vector perpendicular to the cylinder (radius)
perpV = p - repmat(paraM,[1 size(p,2)]).*repmat(cylZ,[size(p,1) 1]);
pOut(:,1) = sqrt(sum(perpV.^2,2));


% Angle of vector around the axis from Z
b0 = [0 0 1];
b0para = cylZ.*b0;
b0proj = b0-b0para;
b0proj = b0proj./sqrt(sum(b0proj.^2));

% IF b0 and cirZ are parallel, any b0proj perpendicular to b0 will do.
if isnan(b0proj)
    b0proj = [1 0 0];
end
b0perp = cross(b0proj,cylZ);
b0perp = b0perp./sqrt(sum(b0perp.^2));

phiParaM = perpV*b0proj';
phiParaV = repmat(phiParaM,[1 size(p,2)]).*repmat(b0proj,[size(p,1) 1]);

phiPerpM = perpV*b0perp';
phiPerpV = repmat(phiPerpM,[1 size(p,2)]).*repmat(b0perp,[size(p,1) 1]);


%phiParaM = sqrt(sum(phiParaV.^2,2));
%phiPerpM = sqrt(sum(phiPerpV.^2,2));

pOut(:,2) = atan(phiPerpM./phiParaM);
q1 = phiPerpM>=0 & phiParaM>=0;
q2 = phiPerpM>=0 & phiParaM<0;
q3 = phiPerpM<0 & phiParaM<0;
q4 = phiPerpM<0 & phiParaM>=0;

pOut(:,2) = pOut(:,2) + pi * (0*q1 + q2 +q3 +2*q4);

%q1 = phiPara>0;

%direction = sqrt(sum((phiPara + repmat(b0proj,[size(p,1) 1])).^2,2))>sqrt(sum(phiPara.^2,2));
%pOut(phiPara>0,2) = 0;

end

