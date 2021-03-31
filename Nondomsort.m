function [Rank,L,f] = Nondomsort(par,rtn)
%Non-Dominated Ranks formation
[N,V] = size(par); F = 2;
Sp = zeros(N,N);
Np = zeros(N,1);
[f,y] = objective(par,rtn,F);
for i=1:N-1
    for j=i+1:N
            comp = f(i,:)-f(j,:);
            n = find(comp<=0);
            if length(n)==2
                Sp(i,j) = 1;
                Np(j) = Np(j) + 1;
            end
            if length(n)<1
                Np(i) = Np(i) + 1;
                Sp(j,i) = 1;
            end
    end
end

Rmax = max(Np);
Rank = zeros(Rmax-1,N);
Nrem = [1:N]; i=1;
while max(Np(Nrem))>=2
    Ri = find(Np(Nrem)==0);
    L(i) = length(Ri);
    Rank(i,1:L(i)) = Nrem(Ri);
    for j=1:L(i)
        P = find(Sp(Nrem(Ri(j)),:)==1);
        Np(P) = Np(P) - 1;
    end
    Nrem = setdiff(Nrem,Nrem(Ri));
    i = i+1;
end
R = i-1;
Rank = Rank(1:R,1:max(L));
end

%for i=1:5
%plot(-f(Rank(i,1:L(i)),1),f(Rank(i,1:L(i)),2),'+')
%hold all
%end
        
 
