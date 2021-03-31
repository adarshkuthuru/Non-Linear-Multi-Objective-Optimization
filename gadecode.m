function [par] = gadecode(pop,nbits,x)
    [N,L] = size(pop); V = size(x,2); xmin = min(x); xmax = max(x);
    par = zeros(N,V);
    for i=1:N
        n = 1;
    for j=1:nbits:L
        B = bi2de(pop(i,j+nbits-1:-1:j));
        par(i,n) = xmin(n) + (xmax(n)-xmin(n))*B/(2^nbits-1);
        n = n + 1;
    end
    end
end