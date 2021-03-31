% function [cost,ind] = bin_GA 
clear all
clc

maxit=200; 
mincost=-9999999; 

x=xlsread('NewData.xlsx','Sheet1','ac2:ag246');
rtn=xlsread('NewData.xlsx','Sheet1','bb2:bf246');
xmin = min(x); xmax = max(x);
npar=size(x,2); 
popsize=size(x,1); % set population size
mutrate=0.1; % set mutation rate
selection=0.5; % fraction of population 
% kept
nbits=10;  % number of bits in each 
% parameter
Nt=nbits*npar;  % total number of bits 
% in a chormosome
keep=floor(selection*popsize); % #population members 
% that survive
pop = [];
popi = zeros(popsize,Nt);
iga=0; % generation counter 
% initialized
for i=1:popsize
    for j=1:npar
        Bi = zeros(1,nbits);
        B = round((x(i,j) - xmin(j))*(2^nbits-1)/(xmax(j)-xmin(j)));
        D = de2bi(B); Bi(1:length(D)) = D;
        popi(i,j*nbits-nbits+1:j*nbits) = fliplr(Bi);
    end
end

par=gadecode(popi,nbits,x); % convert binary to 
% continuous values
[cost,rtnpar]=fitness(par,x,rtn); % calculates population 
% cost using ff
[Rank,L,f] = Nondomsort(par,rtnpar);
[Cd] = Crowding(Rank,f);
for i=1:size(Rank,1)
    pop = [pop;popi(Rank(i,1:L(i)),:)];
end
N = size(pop,1);
if N<popsize
    k = 1; nr = 1;
    while N+k<popsize
        if L(nr)<=popsize-N
        pop(k+N,:) = popi(Rank(nr,k),:);
        k = k+1;
        if k>L(nr)
            nr=nr+1;
        end
        elseif L(nr)>popsize-N
            [C,id] = sort(Cd(nr,1:L(nr)),2,'descend');
            pop(k+N:popsize,:) = popi(Rank(nr,id(1:popsize-N-k+1)),:);
            k=popsize-N+1;
        end
    end   
end
        
% [cost,ind]=sort(cost); % min cost in element 1
% par=par(ind,:);pop=pop(ind,:); % sorts population with 
% lowest cost first
minc(1)=min(cost); % minc contains min of 
% population
meanc(1)=mean(cost); % meanc contains mean 
% of population

while iga<maxit
iga=iga+1; 
M=ceil((popsize-keep)/2); % number of matings
prob=flipud([1:keep]'/sum([1:keep]));% weights 
% chromosomes based 
% upon position in 
% list
odds=[0 cumsum(prob(1:keep))']; % probability
% distribution function
% PROGRAM I: BINARY GENETIC ALGORTHIM 213
pick1=rand(1,M); % mate #1
pick2=rand(1,M); % mate #2
% ma and pa contain the indicies of the chromosomes 
% that will mate
ic=1;
while ic<=M
for id=2:keep+1
if pick1(ic)<=odds(id) && pick1(ic)>odds(id-1)
ma(ic)=id-1;
end % if
if pick2(ic)<=odds(id) && pick2(ic)>odds(id-1)
pa(ic)=id-1;
end % if
end % id
ic=ic+1;
end % while
ix=0:2:keep; % index of mate #1
xp=ceil(rand(1,M)*(Nt-1)); % crossover point
pop(keep+ix,:)=[pop(ma,1:xp) pop(pa,xp+1:Nt)];
% first offspring
pop(keep+ix+1,:)=[pop(pa,1:xp) pop(ma,xp+1:Nt)];
% second offspring
nmut=ceil((popsize-1)*Nt*mutrate); % total number 
% of mutations
mrow=ceil(rand(1,nmut)*(popsize-1))+1; % row to mutate
mcol=ceil(rand(1,nmut)*Nt); % column to mutate
for ii=1:nmut
pop(mrow(ii),mcol(ii))=abs(pop(mrow(ii),mcol(ii))-1);
% toggles bits
end % ii
% The population is re-evaluated for cost
par=gadecode(pop,nbits,x);
% decode
[cost,rtnpar]=fitness(par,x,rtn);
% [cost,ind]=sort(cost);
% par=par(ind,:); pop=pop(ind,:);
minc(iga+1)=min(cost);
meanc(iga+1)=mean(cost);
popn = [];
[Rank,L,f] = Nondomsort(par,rtnpar);
[Cd] = Crowding(Rank,f);
for i=1:size(Rank,1)
    popn = [popn;pop(Rank(i,1:L(i)),:)];
end
N = size(popn,1);
if N<popsize
    k = 1; nr = 1;
    while N+k<popsize
        if L(nr)<=popsize-N
        popn(k+N,:) = pop(Rank(nr,k),:);
        k = k+1;
        if k>L(nr)
            nr=nr+1;
        end
        elseif L(nr)>popsize-N
            [C,id] = sort(Cd(nr,:),2,'descend');
            popn(k+N:popsize,:) = pop(Rank(nr,id(1:popsize-N-k+1)),:);
            k=popsize-N+1;
        end
    end   
end
pop = popn;
if iga>maxit || cost(1)<mincost
break
end
[iga cost(1)]
end 

%plot Ranks
for i=1:size(Rank,1)
    plot(f(Rank(i,1:L(i)),1),f(Rank(i,1:L(i)),2),'o')
    hold all
end
%y=(0.5*f(Rank(1,1:L(1)),1))+(0.5*f(Rank(1,1:L(1)),2));
%[c,i] = min(y);
%x(Rank(1,i),:)