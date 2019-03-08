function [iter,xstar,zstar,opt_basis,nonbasis0,nonbasisu,improving_ray,oporun] = isimplex(A,b,c,u,starting_basis,maxiters,printlevel)
%Impliments a rudimentary simplex method to solve LPs
%   Detailed explanation goes here
basis = starting_basis;
improving_ray = zeros(length(c),1);
count = 1;
tol=10^-8;
tolinf=10^10;
nonbasis0 = setdiff(1:length(c),basis);
nonbasisu = [];
xn=zeros(length(nonbasis0),1);
[rowA,colA] = size(A);

while(count <=maxiters)
    
    
   if mod(count,100) == 0 && printlevel == 2
        fprintf('Current objective value: %d\n',cbB*b + rn*xn)
   end
    
    %assert(min(setdiff(nonbasis0,nonbasisu)==nonbasis0)==1,'error:nonbasis0 and nonbasisu sharing a var');
    %This assert command is only to check to make sure that nonbasis0 and
    %nonbasisu don't contain the same variable. It doesn't work when
    %nonbasis0 is empty but it's just here for testing.
    
    nonbasis = union(nonbasis0,nonbasisu);
    cb = c(basis);
    cn = c(nonbasis);
    B = A(:,basis);
    N = A(:,nonbasis);

    assert(rank(full(B)) == rowA, 'Error, B is singular');
    
    [LB,UB,PB]=lu(B);
    wB=LB\(PB*b);
    xbb = UB\wB;
    
    j=1;
    
    xn=zeros(length(nonbasis),1);
    for(i=1:length(nonbasis))
        if j<=length(nonbasisu) && nonbasis(i) == nonbasisu(j)
            xn(i) = u(nonbasis(i));
            j=j+1;
        end
    end
%    [rowXN,colXN] = size(xn);
%    if rowXN == 1
%        xn = transpose(xn);
%    end
    
    nxn = N*xn;
    wN=LB\(PB*nxn);
    xnb = UB\wN;
    
    
    xb=xbb - xnb;
    xb = tolerance(xb,tol,tolinf,u(basis));
    
    wz = cb/UB;
    vb = wz/LB;
    cbB=vb*PB;
    
    rn = cn - cbB*N;
    rn = tolerance(rn,tol,tolinf);
      
    %Finding the entering variable
    enter1_cost = 0;
    for i=1:length(nonbasis)
        if rn(i) < enter1_cost
            if any(nonbasisu==nonbasis(i)) == 0 
                enter1_cost = rn(i);
                enter1 = i;
            end
        end
    end
    
    
    enter2_cost = 0;
    for i=1:length(nonbasis)
        if rn(i) > enter2_cost
            if any(nonbasis0==nonbasis(i)) == 0
                enter2_cost = rn(i);
                enter2 = i;
            end
        end
    end
    
    if enter1_cost == 0 && enter2_cost == 0
        %optimal solution found
        xstar = zeros(length(c),1);
        xstar(basis)=xb;
	xstar(nonbasisu)=u(nonbasisu);
        zstar = cbB*b + rn*xn;
        opt_basis = basis;
        oporun = 0;
        break
   end
    
    if abs(enter1_cost) >= enter2_cost
        enter_cost = enter1_cost;
        enter_index = enter1;
    else
        enter_cost = enter2_cost;
        enter_index = enter2;
    end
    %Entering variable found   
    %TEST CODE ONLY
    %[enter_cost,enter_index] = min(rn);
    %TEST CODE ONLY
    a=N(:,enter_index);
    wd=LB\(PB*a);
    d = -UB\wd;
    d = tolerance(d,tol,tolinf);
    
    %I think this condition needs to be modified
    if (enter_cost < 0 && min(d) >= 0) || (enter_cost > 0 && max(d) <= 0)
        %unbounded
        xstar = zeros(length(c),1);
        xstar(basis)=xb;
        zstar = -inf;
        opt_basis = basis;
        improving_ray(basis) = d;
        improving_ray(nonbasis(enter_index)) = 1;
        oporun = 1;
        break
    end

    
%    if enter_cost >= 0
%        enter_cost = enter_cost
%    end
    
    if enter_cost < 0
        [leaving_index,bound] = leaving_var(xb,d,u(basis),u(nonbasis(enter_index)),1);
    else
        [leaving_index,bound] = leaving_var(xb,d,u(basis),u(nonbasis(enter_index)),0);
    end
    
    %TEST CODE ONLY
%    if max(d(d<=0)) == 0
%       for n = 1:length(d)
%            if d(n,1) == 0
%                leaving_index = n;
%            end
%        end
%    end
%    xbd = xb./d;
%    leaving=max(xbd(xbd<0));
%    for n = 1:length(xbd)
%        if xbd(n,1) == leaving
%            leaving_index = n;
%            break
%        end
%    end
%    bound = 0;
    %TEST CODE ONLY
    
    
    if leaving_index ~= 0
        ne = nonbasis(enter_index);
        bl = basis(leaving_index);
        basis = union(setdiff(basis,bl),ne);
        nonbasis0 = setdiff(nonbasis0,ne);
        nonbasisu = setdiff(nonbasisu,ne);        
        if bound == 0
            nonbasis0 = union(nonbasis0,bl);
        end
        if bound > 0 && bound < inf
            nonbasisu = union(nonbasisu,bl);
        end        
    else
        if enter_cost < 0
            xn(enter_index) = u(nonbasis(enter_index));
            nonbasis0 = setdiff(nonbasis0,nonbasis(enter_index));
            nonbasisu = union(nonbasisu, nonbasis(enter_index));
        end
        if enter_cost > 0
            xn(enter_index) = 0;
            nonbasis0 = union(nonbasis0,nonbasis(enter_index));
            nonbasisu = setdiff(nonbasisu, nonbasis(enter_index));
        end
    end 

    count=count+1;
end
iter = count;
end