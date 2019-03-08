function [leaving_index_xb,bound] = leaving_var(xb,d,ub,uj,inc_or_dec)
%Determines gives leaving variable index relative to xb. If there is no
%leaving variable, returns 0.
%   Detailed explanation goes here

if inc_or_dec == 1
    if min(xb+uj*d >0) == 1 && min(xb+uj*d<ub) == 1
        leaving_index_xb = 0;
        bound = NaN;
        return;
    end
    t0 = uj;
    tu = uj;
    z = 0;
    u = 0;
    for(i=1:length(xb))
        if d(i) < 0 && -xb(i)/d(i) <= t0
        %if -xb(i)/d(i) <= t0 && -xb(i)/d(i) >= 0
            t0 = -xb(i)/d(i);
            leaving_index0 = i;
            z = z+1;
        end
        if d(i) > 0 && (ub(i)-xb(i))/d(i) <= tu
        %if (ub(i)-xb(i))/d(i) <= tu && (ub(i)-xb(i))/d(i) >= 0
            tu = (ub(i)-xb(i))/d(i);
            leaving_indexu = i;
            u = u+1;
        end
        i = i+1;
    end
    
    t = min(t0,tu);
    %t = t0;
    
    if t == t0 && z > 0
    %if t == t0
        leaving_index_xb = leaving_index0;
        bound = 0;
    elseif t == tu && u > 0
    %elseif t == tu
        leaving_index_xb = leaving_indexu;
        bound = ub(leaving_index_xb);
    else
        disp('oh no why am I here?')
    end
               
end

if inc_or_dec == 0
    if min(xb > 0) == 1 && min(xb < ub) == 1 
            leaving_index_xb = 0;
            bound = NaN;
            return;
    end
    
    t0 = uj;
    tu = uj;
    z = 0;
    u = 0;
    
    for(i = 1:length(xb))
        if d(i) > 0 && xb(i)/d(i) <= t0
            t0 = xb(i)/d(i);
            leaving_index0 = i;
            z = z+1;
        end
        if d(i) < 0 && (xb(i)-ub(i))/d(i) <= tu
                tu = (xb(i)-ub(i))/d(i);
                leaving_indexu = i;
                u = u+1;
        end
        i = i+1;
    end
    
    t = min(t0,tu);
    
    if t == t0 && z > 0
        leaving_index_xb = leaving_index0;
        bound = 0;
    elseif t == tu && u > 0
        leaving_index_xb = leaving_indexu;
        bound = ub(leaving_index_xb);
    else
        disp('oh no why am I here?')
    end
end

end

