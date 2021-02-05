function [param]=inference_g()
%Inference for the metapopulation SEIR model
load Mig_post %load mobility
load pop %load population
load in2020_post %load inner-city mobility data in 2020
load hu % load cityindexs of cities which in Hubei
M=Mig_post(:,:,10:75);
mig=Mig_post(:,:,1:75);
inflow=in2020_post(:,1:66);%.*pop/10e6;
seed=3000;
num_loc=size(M,1);%number of locations
%observation operator: obs=Hx
H=zeros(num_loc,5*num_loc+8);
for i=1:num_loc
    H(i,(i-1)*5+5)=1;
end

load incidence %load observation
incidence=incidence(10:75,:);
%delay day generate
%period one
Td=4;%average reporting delay
a=3.5;%shape parameter of gamma distribution
b=Td/a;%scale parameter of gamma distribution
rnds1=ceil(gamrnd(a,b,1e4,1))-1;%pre-generate gamma random numbers
%period two
Td=2;%average reporting delay
a=3;%shape parameter of gamma distribution
b=Td/a;%scale parameter of gamma distribution
rnds2=ceil(gamrnd(a,b,1e4,1))-1;%pre-generate gamma random numbers
num_times=size(incidence,1);
obs_truth=incidence';
%set OEV
OEV=zeros(num_loc,num_times);
for l=1:num_loc
    for t=1:num_times
        OEV(l,t)=max(4,obs_truth(l,t)^2/4);
    end
end
num_ens=300;%number of ensemble
pop0=pop*ones(1,num_ens);
[x,paramax,paramin]=initialize(pop0,num_ens,seed,M);%get parameter range
num_var=size(x,1);%number of state variables
%IF setting
Iter=10;%number of iterations
num_para=size(paramax,1);%number of parameters
%parameters in different period 
theta=zeros(num_para,Iter+1);%period1: mean parameters at each iteration
theta2=zeros(num_para,Iter+1);%period2: mean parameters at each iteration
theta3=zeros(num_para,Iter+1);%period3: mean parameters at each iteration
theta4=zeros(num_para,Iter+1);%period4: mean parameters at each iteration
theta5=zeros(num_para,Iter+1);%period5: mean parameters at each iteration
para_post=zeros(num_para,num_ens,num_times,Iter);%posterior parameters
sig=zeros(1,Iter);%variance shrinking parameter
alp=0.9;%variance shrinking rate
SIG=(paramax-paramin).^2/4;%initial covariance of parameters
paramax(end-1)=1;
paramin(end-1)=0.05;
SIG2=(paramax-paramin).^2/4;%initial covariance of parameters
lambda=1.1;%inflation parameter to aviod divergence within each iteration
%start iteration for Iter round
param=zeros(num_para,num_ens,num_times,Iter);
for n=1:Iter
    sig(n)=alp^(n-1);
    %generate new ensemble members using multivariate normal distribution
    Sigma=diag(sig(n)^2*SIG);
    Sigma2=diag(sig(n)^2*SIG2);
    if (n==1)
        %first guess of state space
        [x,~,~]=initialize(pop0,num_ens,seed,M);
        para=x(end-num_para+1:end,:);
        theta(:,1)=mean(para,2);%mean parameter
        theta2(:,1)=mean(para,2);%mean parameter
        theta3(:,1)=mean(para,2);%mean parameter
        theta4(:,1)=mean(para,2);%mean parameter
        theta5(:,1)=mean(para,2);%mean parameter
    else
        [x,~,~]=initialize(pop0,num_ens,seed,M);
        para=mvnrnd(theta(:,n)',Sigma,num_ens)';%generate parameters
        x(end-num_para+1:end,:)=para;
    end
    %correct lower/upper bounds of the parameter
    x=checkbound_ini(x,pop0);
    %Begin looping through observations
    x_prior=zeros(num_var,num_ens,num_times);%prior
    x_post=zeros(num_var,num_ens,num_times);%posterior
    pop=pop0;
    obs_temp=zeros(num_loc,num_ens,num_times);%records of reported cases
    for t=1:num_times
        [n,t]
        if n>1
        if t==15
            para=mvnrnd(theta2(:,n)',Sigma,num_ens)';%generate parameters
            x(end-num_para+1:end,:)=para;
        elseif t==25
            para=mvnrnd(theta3(:,n)',Sigma2,num_ens)';%generate parameters
            x(end-num_para+1:end,:)=para;
            %elseif t==30
            %    para=mvnrnd(theta4(:,n)',Sigma,num_ens)';%generate parameters
            %    x(end-num_para+1:end,:)=para;
        elseif t==35
            para=mvnrnd(theta4(:,n)',Sigma2,num_ens)';%generate parameters
            x(end-num_para+1:end,:)=para;
        elseif t==45
            para=mvnrnd(theta5(:,n)',Sigma2,num_ens)';%generate parameters
            x(end-num_para+1:end,:)=para;
        end
        else
            if t==35
                para=mvnrnd(theta(:,n)',Sigma,num_ens)';%generate parameters
                x(end,:)=para(5,:);
            end
        end
        %inflation
        x=mean(x,2)*ones(1,num_ens)+lambda*(x-mean(x,2)*ones(1,num_ens));
        if t<34 
           x = checkbound2(x,pop);
        else
           x = checkbound(x,pop);
        end
        %integrate forward
        [x,pop]=SEIR(x,M,inflow,pop,t,pop0);
        obs_cnt=H*x;%new infection
        %add reporting delay
        for k=1:num_ens
            for l=1:num_loc
                if obs_cnt(l,k)>0
                   if (t==34 && ismember(l,hu)) %12feb2020 no delay documented infected
                       rnds=0;
                   elseif t<23
                       rnds=rnds1;
                   else
                       rnds=rnds2;
                   end
                   rnd=datasample(rnds,obs_cnt(l,k));
                    for h=1:length(rnd)
                        if (t+rnd(h)<=num_times)
                            obs_temp(l,k,t+rnd(h))=obs_temp(l,k,t+rnd(h))+1;
                        end
                    end
                    if (t==34 && ismember(l,hu))
                       obs_temp(l,k,34)=squeeze(sum(obs_temp(l,k,34:end),3));
                       obs_temp(1,k,35:end)=0;
                    end
                end
            end
        end
        obs_ens=obs_temp(:,:,t);%observation at t
        x_prior(:,:,t)=x;%set prior
        %loop through local observations
        for l=1:num_loc
            %Get the variance of the ensemble
            obs_var = OEV(l,t);
            prior_var = var(obs_ens(l,:));
            post_var = prior_var*obs_var/(prior_var+obs_var);
            if prior_var==0 %if degenerate
                post_var=1e-3;
                prior_var=1e-3;
            end
            prior_mean = mean(obs_ens(l,:));
            post_mean = post_var*(prior_mean/prior_var + obs_truth(l,t)/obs_var);
            %%%% Compute alpha and adjust distribution to conform to posterior moments
            alpha = (obs_var/(obs_var+prior_var)).^0.5;
            dy = post_mean + alpha*(obs_ens(l,:)-prior_mean)-obs_ens(l,:);
            %Loop over each state variable (connected to location l)
            rr=zeros(1,num_var);
            
            jant=t+9;
            neighbors=union(find(sum(mig(l,:,jant-7:jant))>0),find(sum(mig(:,l,jant-7:jant))>0));
            neighbors=[neighbors;l];%add location l
            for i=1:length(neighbors)
                idx=neighbors(i);
                for j=1:5
                   A=cov(x((idx-1)*5+j,:),obs_ens(l,:));
                   rr((idx-1)*5+j)=A(2,1)/prior_var;
                end
            end
            for i=num_loc*5+1:num_loc*5+7
                A=cov(x(i,:),obs_ens(l,:));
                rr(i)=A(2,1)/prior_var;
            end
            if ismember(l,hu)
               for i=num_loc*5+8
                   A=cov(x(i,:),obs_ens(l,:));
                   rr(i)=A(2,1)/prior_var;
               end
            end
            %Get the adjusted variable
            dx=rr'*dy;
            x=x+dx;
            %Corrections to DA produced aphysicalities
            if t<34
                x = checkbound2(x,pop);
            else
                x = checkbound(x,pop);
            end
        end
        x_post(:,:,t)=x;
        para_post(:,:,t,n)=x(end-num_para+1:end,:);
    end
    para=x_post(end-num_para+1:end,:,1:num_times);
    temp=squeeze(median(para,2));%average over ensemble members
    theta(:,n+1)=mean(temp(:,1:14),2);%average over time
    theta2(:,n+1)=mean(temp(:,15:24),2);%average over time
    theta3(:,n+1)=mean(temp(:,27:34),2);%average over time
    theta4(:,n+1)=mean(temp(:,35:44),2);%average over time
    theta5(:,n+1)=mean(temp(:,45:end),2);%average over time
    param(:,:,:,n)=para;
end

save('para2','param');
end

function x = checkbound_ini(x,pop)
%S,E,Is,Ia,obs,...,beta,mu,theta,Z,alpha,D
betalow=0.22;betaup=0.5;%transmission rate
mulow=0.05;muup=1.0;%relative transmissibility
Zlow=2;Zup=5;%latency period
alphalow=0.05;alphaup=1;%reporting rate
alphalow2=0.01;alphaup2=0.1;%reporting rate
D2low=3;D2up=5;%infectious period
Dlow=1;Dup=3;%infectious period
xmin=[mulow;betalow;mulow;Zlow;alphalow;Dlow;D2low;alphalow2];
xmax=[muup;betaup;muup;Zup;alphaup;Dup;D2up;alphaup2];
num_loc=size(pop,1);
for i=1:num_loc
    %S
    x((i-1)*5+1,x((i-1)*5+1,:)<0)=0;
    x((i-1)*5+1,x((i-1)*5+1,:)>pop(i,:))=pop(i,x((i-1)*5+1,:)>pop(i,:));
    %E
    x((i-1)*5+2,x((i-1)*5+2,:)<0)=0;
    %Ir
    x((i-1)*5+3,x((i-1)*5+3,:)<0)=0;
    %Iu
    x((i-1)*5+4,x((i-1)*5+4,:)<0)=0;
    %obs
    x((i-1)*5+5,x((i-1)*5+5,:)<0)=0;
end
for i=1:8
    temp=x(end-8+i,:);
    index=(temp<xmin(i))|(temp>xmax(i));
    index_out=find(index>0);
    index_in=find(index==0);
    %redistribute out bound ensemble members
    x(end-8+i,index_out)=datasample(x(end-8+i,index_in),length(index_out));
end
end

function x = checkbound(x,pop)
%S,E,Is,Ia,obs,...,beta,mu,theta,Z,alpha,D
betalow=0.22;betaup=0.5;%transmission rate
mulow=0.05;muup=1.0;%relative transmissibility
Zlow=2;Zup=5;%latency period
alphalow=0.05;alphaup=1.0;%reporting rate
D2low=3;D2up=5;%infectious period
Dlow=1;Dup=3;%infectious period
xmin=[mulow;betalow;mulow;Zlow;alphalow;Dlow;D2low;alphalow];
xmax=[muup;betaup;muup;Zup;alphaup;Dup;D2up;alphaup];
num_loc=size(pop,1);
for i=1:num_loc
    %S
    x((i-1)*5+1,x((i-1)*5+1,:)<0)=0;
    x((i-1)*5+1,x((i-1)*5+1,:)>pop(i,:))=pop(i,x((i-1)*5+1,:)>pop(i,:));
    %E
    x((i-1)*5+2,x((i-1)*5+2,:)<0)=0;
    %Ir
    x((i-1)*5+3,x((i-1)*5+3,:)<0)=0;
    %Iu
    x((i-1)*5+4,x((i-1)*5+4,:)<0)=0;
    %obs
    x((i-1)*5+5,x((i-1)*5+5,:)<0)=0;
end
for i=1:8
    x(end-8+i,x(end-8+i,:)<xmin(i))=xmin(i)*(1+0.1*rand(sum(x(end-8+i,:)<xmin(i)),1));
    x(end-8+i,x(end-8+i,:)>xmax(i))=xmax(i)*(1-0.1*rand(sum(x(end-8+i,:)>xmax(i)),1));
end
end 

function x = checkbound2(x,pop)
%S,E,Is,Ia,obs,...,beta,mu,theta,Z,alpha,D
betalow=0.22;betaup=0.5;%transmission rate
mulow=0.05;muup=1.0;%relative transmissibility
Zlow=2;Zup=5;%latency period
alphalow=0.05;alphaup=1.0;%reporting rate
alphalow2=0.01;alphaup2=0.1;%reporting rate
D2low=3;D2up=5;%infectious period
Dlow=1;Dup=3;%infectious period
xmin=[mulow;betalow;mulow;Zlow;alphalow;Dlow;D2low;alphalow2];
xmax=[muup;betaup;muup;Zup;alphaup;Dup;D2up;alphaup2];
num_loc=size(pop,1);
for i=1:num_loc
    %S
    x((i-1)*5+1,x((i-1)*5+1,:)<0)=0;
    x((i-1)*5+1,x((i-1)*5+1,:)>pop(i,:))=pop(i,x((i-1)*5+1,:)>pop(i,:));
    %E
    x((i-1)*5+2,x((i-1)*5+2,:)<0)=0;
    %Ir
    x((i-1)*5+3,x((i-1)*5+3,:)<0)=0;
    %Iu
    x((i-1)*5+4,x((i-1)*5+4,:)<0)=0;
    %obs
    x((i-1)*5+5,x((i-1)*5+5,:)<0)=0;
end
for i=1:8
    x(end-8+i,x(end-8+i,:)<xmin(i))=xmin(i)*(1+0.1*rand(sum(x(end-8+i,:)<xmin(i)),1));
    x(end-8+i,x(end-8+i,:)>xmax(i))=xmax(i)*(1-0.1*rand(sum(x(end-8+i,:)>xmax(i)),1));
end
end
