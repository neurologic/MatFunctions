function f = ImpulseFunction(xtime,s)
%s is spike times (1,s) vector

f = zeros(1,size(xtime,2));


for i = 1 : size(s,2)
    
    t = xtime - s(i);
    % Our default value is 0
    inds = find(t ==0);
    
    if size(inds,2)>1
        f(inds(1)) = 1;
        
    else
        f(inds) = 1;
    end
    
    
    
end
