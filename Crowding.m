function [Cd] = Crowding(Rank,f)
    [R,L] = size(Rank);
    m = size(f,2); Cd = zeros(R,L);
    Cd(:,1) = 100;
    for i=1:R
        Ndnum = find(Rank(i,:)>0);
        Ndset = Rank(i,Ndnum);
        for k=2:length(Ndset)-1
            Cd(i,length(Ndset)) = 100;
        for j=1:m
            [F1,I] = sort(f(Ndset,j));
            Cd(i,k) = Cd(i,k) + abs(f(I(k+1),j)-f(I(k-1),j))/abs(F1(end)-F1(1));
        end
        end
    end          
end
