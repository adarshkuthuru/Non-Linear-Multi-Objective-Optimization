function [f,y]=objective(x,rtn,k)
[N,V] = size(x);
f = zeros(N,k);
% Objective function one - Returns
for j=1:N
sum1=0;
for i = 1:V
    sum1 = sum1+(x(j,i)*rtn(j,i));
end
% Decision variables are used to form the objective function.
f(j,k-1) = -sum1;

% Objective function two - CVaR
%Numerator of CVaR
ret1=0;ret2=0;ret3=0;ret4=0;ret5=0;n=floor(0.05*N);
for i=1:n
    ret1 = ret1 + rtn(i,1);
    ret2 = ret2 + rtn(i,2);
    ret3 = ret3 + rtn(i,3);
    ret4 = ret4 + rtn(i,4);
    ret5 = ret5 + rtn(i,5);
end
r(1) = ret1/n;r(2) = ret2/n;r(3) = ret3/n;r(4) = ret4/n;r(5) = ret5/n;

sum2=0;
for i = 1:V
    sum2 = sum2+(x(j,i)*r(i));
end
% Decision variables are used to form the objective function.
f(j,k) = abs(sum2);

y(j)=(0.5*f(j,k-1))+(0.5*f(j,k));
end
end


