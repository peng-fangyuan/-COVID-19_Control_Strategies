function [x,paramax,paramin]=initialize(pop,num_ens,seed,M)
load Mig %load mobility
%Initialize the metapopulation SEIRS model
num_loc=size(pop,1);
% num_var=5*num_loc+6;
% S,E,Is,Ia,obs,...,beta,mu,theta,Z,alpha,D
% prior range
Slow=1.0;Sup=1.0;%susceptible fraction
Elow=0;Eup=0;%exposed
Irlow=0;Irup=0;%documented infection
Iulow=0;Iuup=0;%undocumented infection
obslow=0;obsup=0;%reported case
betalow=0.22;betaup=0.5;%transmission rate
mulow=0.05;muup=1.0;%relative transmissibility
Zlow=2;Zup=5;%latency period
alphalow=0.05;alphaup=1.0;%reporting rate
alphalow2=0.01;alphaup2=0.1;%reporting rate
Dlow=1;Dup=3;%infectious period
D2low=2;D2up=5;
%range of model state including variables and parameters
xmin=[];
xmax=[];
for i=1:num_loc
    xmin=[xmin;Slow*pop(i);Elow*pop(i);Irlow*pop(i);Iulow*pop(i);obslow];
    xmax=[xmax;Sup*pop(i);Eup*pop(i);Irup*pop(i);Iuup*pop(i);obsup];
end
xmin=[xmin;mulow;betalow;mulow;Zlow;alphalow;Dlow;D2low;alphalow2];
xmax=[xmax;muup;betaup;muup;Zup;alphaup;Dup;D2up;alphaup2];
paramax=xmax(end-7:end);
paramin=xmin(end-7:end);

%seeding in Wuhan
%Wuhan 
seedid=292;
%E
xmin((seedid-1)*5+2)=0;xmax((seedid-1)*5+2)=seed;
%Is
xmin((seedid-1)*5+3)=0;xmax((seedid-1)*5+3)=0;
%Ia
xmin((seedid-1)*5+4)=0;xmax((seedid-1)*5+4)=seed;
%Latin Hypercubic Sampling
x=lhsu(xmin,xmax,num_ens);
x=x';
for i=1:num_loc
    x((i-1)*5+1:(i-1)*5+4,:)=round(x((i-1)*5+1:(i-1)*5+4,:));
end
%seeding in other cities
C=M(:,292,1);%first day
%C=sum(mig(:,292,5:10),3);%first day
for i=1:num_loc
    if i~=seedid
        %E
        Ewuhan=x((seedid-1)*5+2,:);
        x((i-1)*5+2,:)=round(C(i)*6.5*Ewuhan/pop(seedid));
        %Ia
        Iawuhan=x((seedid-1)*5+4,:);
        x((i-1)*5+4,:)=round(C(i)*6.5*Iawuhan/pop(seedid));
    end
end

function s=lhsu(xmin,xmax,nsample)
nvar=length(xmin);
ran=rand(nsample,nvar);
s=zeros(nsample,nvar);
for j=1: nvar
   idx=randperm(nsample);
   P =(idx'-ran(:,j))/nsample;
   s(:,j) = xmin(j) + P.* (xmax(j)-xmin(j));
end