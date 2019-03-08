function [m] = tolerance(m,tol,tolinf,u)
%Corrects for tolerances
%   Detailed explanation goes here
if nargin <= 3
    u = inf(size(m));
end
for i = 1:length(m)
    if m(i) >= -tol && m(i) <= tol
        m(i) = 0;
    end
    if m(i) >= u(i)-tol && m(i) <= u(i)+tol
        m(i) = u(i);
    end
    if m(i) >= tolinf
        m(i) = inf;
    end

end

