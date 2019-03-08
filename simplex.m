function [xstar,zstar,opt_basis,nonbasis0,nonbasisu,improving_ray] = simplex(A,b,c,u,maxiters,printlevel)
%Impliments a less rudimentary simplex method to solve LPs
%   Detailed explanation goes here
if nargin <= 4
    maxiters = 1000;
    printlevel = 0;
end
A=sparse(A);
[A,b,p,rank]=rreqns(A,b);
[A,b,c,u,r,s]=scalelp(A,b,c,u);

[rowc,colc]=size(c);
if(rowc ~= 1)
    c = transpose(c);
end

%Phase 1
[rowA,colA] = size(A);
c1 = horzcat(zeros(1,length(c)),ones(1,rowA));
A1 = horzcat(A,zeros(rowA));
    for n = 1:rowA
        if b(n) >= 0
            A1(n,colA+n) = 1;
        else
            A1(n,colA+n) = -1;
        end
    end
u1=vertcat(u,inf(rowA,1));
[iter1,x,w,p1basis] = isimplex(A1,b,c1,u1,[(length(c)+1):length(c1)],maxiters,printlevel);
    if printlevel >= 1
        disp('Phase 1');
        fprintf('iterations: %d\n',iter1);
        fprintf('optimal value: %d\n',w);
    end
    if w ~= 0
        %infeasible
        xstar=NaN*ones(length(c),1);
        zstar = NaN;
        opt_basis=NaN*ones(1,rowA);
        nonbasis0=NaN*ones(1,length(c)-rowA);
        nonbasisu=[];
        improving_ray=NaN*ones(length(c),1);
        if printlevel >= 1
            disp('LP is infeasible');
        end
        return
    end
%    for n = (colA+1):length(x)
 %       assert(x(n) ~= 0, 'Error: Something weird happened in Phase 1')
  %  end
%end Phase 1

%Phase 2
%disp('phase 2')

[rowc,colc]=size(c);
if(rowc ~= 1)
    c = transpose(c);
end

[iter2,xstar,zstar,opt_basis,nonbasis0,nonbasisu,improving_ray,oporun] = isimplex(A,b,c,u,p1basis,maxiters,printlevel);
xstar=unscalesoln(xstar,r,s);
if printlevel >= 1
    disp('Phase 2')
    fprintf('iterations: %d\n',iter2);
    fprintf('optimal value: %d\n',zstar);
    if oporun == 0
        disp('LP is optimal');
    else
        disp('LP is unbounded');
    end
end
%End Phase 2
end

