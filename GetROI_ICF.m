function [ cylROI ] = GetROI_ICF( cylMap )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

cylNums = unique(cylMap(cylMap>0));

cylROI(1:numel(cylNums)) = struct('Id',0,'Bounds',zeros(3,2),'Dir',0);

for i=1:numel(cylNums)
    cylROI(i).Id = cylNums(i);
    
    [x,y,z] = ind2sub(size(cylMap),find(cylMap==cylNums(i)));
    cylROI(i).Bounds(1,:) = [min(x) max(x)];
    cylROI(i).Bounds(2,:) = [min(y) max(y)];
    cylROI(i).Bounds(3,:) = [min(z) max(z)];
    
    [~,cylROI(i).Dir] = max(diff(cylROI(i).Bounds,1,2));
    
    % Hard-code for experiments
    if cylNums(i)==30
        cylROI(i).Dir=3;
    end
    
    cylROI(i).Bounds = cylROI(i).Bounds + ([1 2 3]'~=cylROI(i).Dir)*[-3 3];
end

end

