function [ t ] = TSinT_ICF( y )
%TSINT_VP Summary of this function goes here
%   Detailed explanation goes here

nIterations = 20;

t1 = 0;
t2 = pi;
y1 = TSinTFnc_ICF(t1);
y2 = TSinTFnc_ICF(t2);

for i=1:nIterations
    t = .5*(t1+t2);
    f = TSinTFnc_ICF(t);
    
    if f<y
        y1 = f;
        t1 = t;
    else
        y2 = f;
        t2 = t;
    end
end

t = t1 + (t2-t1)*(y2-y)/(y2-y1);

end

