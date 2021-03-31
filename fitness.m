function [cost,rtnpar] = fitness(par,x,rtn)
    [N,V] = size(par);
    rtnpar = zeros(N,V);
    for j=1:N
    for i=1:V
        Hi = find(x(:,i)>=par(j,i));
        Lo = find(x(:,i)<=par(j,i));
        [H,iH] = sort(x(Hi,i));
        [L,iL] = sort(x(Lo,i),1,'descend');
        rtnpar(j,i) = (rtn(iH(1),i)+rtn(iL(1),i))/2;
    end
    end
    for j=1:N
        sum1=0;
    for i = 1:V
        sum1 = sum1+(par(j,i)*rtnpar(j,i));
    end
    % Decision variables are used to form the objective function.
    y1 = -sum1;
    ret1=0;ret2=0;ret3=0;ret4=0;ret5=0;n=floor(0.05*N);
    for i=1:n
        ret1 = ret1 + rtnpar(i,1);
        ret2 = ret2 + rtnpar(i,2);
        ret3 = ret3 + rtnpar(i,3);
        ret4 = ret4 + rtnpar(i,4);
        ret5 = ret5 + rtnpar(i,5);
    end
    r(1) = ret1/n;r(2) = ret2/n;r(3) = ret3/n;r(4) = ret4/n;r(5) = ret5/n;

    sum2=0;
    for i = 1:V
        sum2 = sum2+(par(j,i)*r(i));
    end
    % Decision variables are used to form the objective function.
    y2 = sum2;
    cost(j) = 0.5*(y1+y2);
    end
end