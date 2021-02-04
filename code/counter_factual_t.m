function counter_factual_t()
load pa_t2

infectall_30LA=zeros(365,188,200);
infectall_50LA=zeros(365,188,200);
infectall_70LA=zeros(365,188,200);
infectall_90LA=zeros(365,188,200);

infectall_30L=zeros(365,188,200);
infectall_50L=zeros(365,188,200);
infectall_70L=zeros(365,188,200);
infectall_90L=zeros(365,188,200);
%infectall_90L=zeros(365,188,200);


parfor i=1:200
    pa=pa_t(:,i);
    %%
    %Only testing and treatment
    [~,infect_30]=optim(pa,0,0,188,0.3); % testing rate 30%
    [~,infect_50]=optim(pa,0,0,188,0.5); % testing rate 50%  
    [~,infect_70]=optim(pa,1,0,188,0.7);  %testing rate 70%  
    [~,infect_90]=optim(pa,1,0,188,0.9); % testing rate 90%
    [~,infect_90]=optim(pa,0,0,188,0.9); % testing rate 90%  
    %%
    %T&T with Lockdown
    [~,infect_30L]=optim(pa,1,0,188,0.3); % testing rate 30%
    [~,infect_50L]=optim(pa,1,0,188,0.5); % testing rate 50%
    [~,infect_70L]=optim(pa,1,0,188,0.7); % testing rate 70%
    [~,infect_90L]=optim(pa,1,0,188,0.9); % testing rate 90%
    %%
    % All policy
    [~,infect_30LA]=optim(pa,1,1,188,0.3); % testing rate 50%
    [~,infect_50LA]=optim(pa,1,1,188,0.5); % testing rate 50%
    [~,infect_70LA]=optim(pa,1,1,188,0.7); % testing rate 70%
    [~,infect_90LA]=optim(pa,1,1,188,0.9); % testing rate 90%

    
end


end

%%
function [obs,infect]=optim(param,intype,migtype,max_t,rate)
%load M
load in2020_post
load in2019_post
load Mig2019_post
load Mig_post
if intype==1
    inflow=in2020_post;% innercity activities in 2020
else
    inflow=in2019_post; % innercity activities in 2019
end
if migtype==1
    M=Mig_post; %intercity activities in 2020
else
    M=Mig2019_post; % intercity activities in 2019
end
num_times=max_t;
[obs,infect]=SEIR_optim(param,inflow,M,num_times,rate);
obs=obs';
%save('obs','obs')
end
function [obs,infect]=SEIR_optim(param,inflow,M,num_times,rate)
load pop
seed=3000;
x_ini=initialize_b(pop,M,seed);
load hu
Td=4;%average reporting delay
a=3.5;%shape parameter of gamma distribution
b=Td/a;%scale parameter of gamma distribution
rnds1=ceil(gamrnd(a,b,1e4,1))-1;
Td=2;%average reporting delay
a=3;%shape parameter of gamma distribution
b=Td/a;%scale parameter of gamma distribution
rnds2=ceil(gamrnd(a,b,1e4,1))-1;%pre-generate gamma random numbers
num_loc=size(M,1);%number of locations
H=zeros(num_loc,6*num_loc);
for i=1:num_loc
    H(i,(i-1)*6+5)=1;
end
H2=zeros(num_loc,6*num_loc);
for i=1:num_loc
    H2(i,(i-1)*6+6)=1;
end
obs=zeros(num_loc,num_times);%records of reported cases
infect=zeros(num_loc,num_times);%everyday new cases
x=x_ini;
pop0=pop;
for t=1:num_times
    [x,pop,M]=SEIR_a2(x,M,inflow,pop,t,pop0,param,rate,infect);
    obs_cnt=H*x;%new infection
    infect(:,t)=H2*x;
    for l=1:num_loc       
        if obs_cnt(l)>0
           if (t==34  && ismember(l,hu))
                       rnds=rnds2;
           elseif t<23
                       rnds=rnds1;
           else
                       rnds=rnds2;
           end
           rnd=datasample(rnds,obs_cnt(l));
           for h=1:length(rnd)
               if (t+rnd(h)<=num_times)
                   obs(l,t+rnd(h))=obs(l,t+rnd(h))+1;
               end               
           end
           if (t==34 && ismember(l,hu))
               obs(l,34)=sum(obs(l,34:end),2);
           %    obs(l,35:end)=0;
           end
        end        
    end
end

end

function [x,pop,M]=SEIR_a2(x,M,innerflow,pop,ts,pop0,param,rate,infect)
%the metapopulation SEIR model
load hu

dt=1;
tmstep=1;
%integrate forward for one day
num_loc=size(pop,1);

%S,E,Is,Ia,obs,...,beta,mu,theta,Z,alpha,D
Sidx=(1:6:6*num_loc)';
Eidx=(2:6:6*num_loc)';
Isidx=(3:6:6*num_loc)';
Iaidx=(4:6:6*num_loc)';
obsidx=(5:6:6*num_loc)';
infidx=(6:6:6*num_loc)';

S=zeros(num_loc,tmstep+1);
E=zeros(num_loc,tmstep+1);
Is=zeros(num_loc,tmstep+1);
Ia=zeros(num_loc,tmstep+1);
Incidence=zeros(num_loc,tmstep+1);
obs=zeros(num_loc,1/dt);
infection=zeros(num_loc,tmstep+1);
%initialize S,E,Is,Ia and parameters
S(:,1)=x(Sidx);
E(:,1)=x(Eidx);
Is(:,1)=x(Isidx);
Ia(:,1)=x(Iaidx);



para=param(1:8);
%with treatment + gradually
if ts<15
    %para=param(1:10);
    
    if rate<param(5)
        para(5)=rate;
    else
        para(5)=param(5);
    end
    para(8)=param(8)*param(5);
elseif ts<25
    %para=param(11:20);
    para(1)=param(9);
    para(6)=param(14);
    if rate<param(5)
        para(5)=rate;
    else
        para(5)=param(5);
    end
    para(8)=param(16)*param(13);
    para(3)=param(11);
elseif ts<34   
    para(1)=param(17);
    para(6)=param(22);
    para(5)=rate;
    para(8)=param(24)*param(21);
    para(3)=param(11);
elseif ts<45
    para(5)=rate;
    para(8)=rate;
    para(1)=param(25);
    para(6)=param(30);
    para(3)=param(11);
else
    para(5)=rate;
    para(8)=rate;
    para(1)=param(33);
    para(6)=param(38);
    para(3)=param(11);
end

mu=ones(num_loc,1)*param(1);
mu2=ones(num_loc,1)*param(1);
theta=6.5;

%alpha and return rate

if ts<180
    Migrant=squeeze(M(:,:,ts));
    inflow=innerflow(:,ts);
else
    Migrant=squeeze(M(:,:,180));
    inflow=innerflow(:,180);
end

tcnt=0;
mu=ones(num_loc,1)*para(3);
mu2=ones(num_loc,1)*para(3);
Z=para(4);

D=ones(num_loc,1)*para(7);
D2=ones(num_loc,1)*para(7);


alpha=ones(num_loc,1)*para(5);
alpha2=ones(num_loc,1)*para(5);

betas=inflow*para(1)*para(2);
betaa=inflow*para(2);

alpha(292)=para(8);
alpha2(292)=para(8);



for t=ts+dt:dt:ts+tmstep
    tcnt=tcnt+1;
    dt1=dt;
    %first step
    ESenter=dt1*(ones(num_loc,1)*theta).*(Migrant*(S(:,tcnt)./(pop-Is(:,tcnt)-theta*(1-mu).*Ia(:,tcnt))));
    ESleft=min(dt1*(ones(num_loc,1)*theta).*(S(:,tcnt)./(pop-Is(:,tcnt)-theta*(1-mu).*Ia(:,tcnt))).*(sum(Migrant)'),dt1*S(:,tcnt));
    EEenter=dt1*(ones(num_loc,1)*theta.*mu2).*(Migrant*(E(:,tcnt)./(pop-Is(:,tcnt)-theta*(1-mu).*Ia(:,tcnt))));
    EEleft=min(dt1*(ones(num_loc,1)*theta.*mu2).*(E(:,tcnt)./(pop-Is(:,tcnt)-theta*(1-mu).*Ia(:,tcnt))).*(sum(Migrant)'),dt1*E(:,tcnt));
    EIaenter=dt1*(ones(num_loc,1)*theta.*mu).*(Migrant*(Ia(:,tcnt)./(pop-Is(:,tcnt)-theta*(1-mu).*Ia(:,tcnt))));
    EIaleft=min(dt1*(ones(num_loc,1)*theta.*mu).*(Ia(:,tcnt)./(pop-Is(:,tcnt)-theta*(1-mu).*Ia(:,tcnt))).*(sum(Migrant)'),dt1*Ia(:,tcnt));
    
    Eexps=dt1*(betas).*S(:,tcnt).*Is(:,tcnt)./pop;
    Eexpa=dt1*(betaa).*S(:,tcnt).*Ia(:,tcnt)./pop;
    Einfs=dt1*(alpha).*E(:,tcnt)./(ones(num_loc,1)*Z);
    Einfa=dt1*((1-alpha)).*E(:,tcnt)./(ones(num_loc,1)*Z);
    Erecs=dt1*Is(:,tcnt)./(D);
    Etoia=dt1*(alpha2).*Ia(:,tcnt);
    Ereca=dt1*Ia(:,tcnt)./(D2);
    
    ESenter=max(ESenter,0);ESleft=max(ESleft,0);
    EEenter=max(EEenter,0);EEleft=max(EEleft,0);
    EIaenter=max(EIaenter,0);EIaleft=max(EIaleft,0);
    Eexps=max(Eexps,0);Eexpa=max(Eexpa,0);
    Einfs=max(Einfs,0);Einfa=max(Einfa,0);
    Erecs=max(Erecs,0);Ereca=max(Ereca,0);
    Etoia=max(Etoia,0);
    
    %%%%%%%%%%stochastic version
    ESenter=poissrnd(ESenter);ESleft=poissrnd(ESleft);
    EEenter=poissrnd(EEenter);EEleft=poissrnd(EEleft);
    EIaenter=poissrnd(EIaenter);EIaleft=poissrnd(EIaleft);
    Eexps=poissrnd(Eexps);
    Eexpa=poissrnd(Eexpa);
    Einfs=poissrnd(Einfs);
    Einfa=poissrnd(Einfa);
    Erecs=poissrnd(Erecs);
    Ereca=poissrnd(Ereca);
    Etoia=poissrnd(Etoia);

    sk1=-Eexps-Eexpa+ESenter-ESleft;
    ek1=Eexps+Eexpa-Einfs-Einfa+EEenter-EEleft;
    isk1=Einfs-Erecs+Etoia;
    iak1=Einfa+EIaenter-EIaleft-Etoia-Ereca;
    ik1i=Einfs+Etoia;
    if1=Einfs+Einfa;
    %second step
    Ts1=S(:,tcnt)+sk1/2;
    Te1=E(:,tcnt)+ek1/2;
    Tis1=Is(:,tcnt)+isk1/2;
    Tia1=Ia(:,tcnt)+iak1/2;
    
    ESenter=dt1*(ones(num_loc,1)*theta).*(Migrant*(Ts1./(pop-Tis1-theta*(1-mu).*Tia1)));
    ESleft=min(dt1*(ones(num_loc,1)*theta).*(Ts1./(pop-Tis1-theta*(1-mu).*Tia1)).*(sum(Migrant)'),dt1*Ts1);
    EEenter=dt1*(ones(num_loc,1)*theta.*mu2).*(Migrant*(Te1./(pop-Tis1-theta*(1-mu).*Tia1)));
    EEleft=min(dt1*(ones(num_loc,1)*theta.*mu2).*(Te1./(pop-Tis1-theta*(1-mu).*Tia1)).*(sum(Migrant)'),dt1*Te1);
    EIaenter=dt1*(ones(num_loc,1)*theta.*mu).*(Migrant*(Tia1./(pop-Tis1-theta*(1-mu).*Tia1)));
    EIaleft=min(dt1*(ones(num_loc,1)*theta.*mu).*(Tia1./(pop-Tis1-theta*(1-mu).*Tia1)).*(sum(Migrant)'),dt1*Tia1);
    
    Eexps=dt1*(betas).*Ts1.*Tis1./pop;
    Eexpa=dt1*(betaa).*Ts1.*Tia1./pop;
    Einfs=dt1*(alpha).*Te1./(ones(num_loc,1)*Z);
    Einfa=dt1*((1-alpha)).*Te1./(ones(num_loc,1)*Z);
    Erecs=dt1*Tis1./(D);
    Etoia=dt1*(alpha2).*Tia1;
    Ereca=dt1*Tia1./(D2);
    
    ESenter=max(ESenter,0);ESleft=max(ESleft,0);
    EEenter=max(EEenter,0);EEleft=max(EEleft,0);
    EIaenter=max(EIaenter,0);EIaleft=max(EIaleft,0);
    Eexps=max(Eexps,0);Eexpa=max(Eexpa,0);
    Einfs=max(Einfs,0);Einfa=max(Einfa,0);
    Erecs=max(Erecs,0);Ereca=max(Ereca,0);
    Etoia=max(Etoia,0);
    
    %%%%%%%%%%stochastic version
    ESenter=poissrnd(ESenter);ESleft=poissrnd(ESleft);
    EEenter=poissrnd(EEenter);EEleft=poissrnd(EEleft);
    EIaenter=poissrnd(EIaenter);EIaleft=poissrnd(EIaleft);
    Eexps=poissrnd(Eexps);
    Eexpa=poissrnd(Eexpa);
    Einfs=poissrnd(Einfs);
    Einfa=poissrnd(Einfa);
    Erecs=poissrnd(Erecs);
    Ereca=poissrnd(Ereca);
    Etoia=poissrnd(Etoia);

    sk2=-Eexps-Eexpa+ESenter-ESleft;
    ek2=Eexps+Eexpa-Einfs-Einfa+EEenter-EEleft;
    isk2=Einfs-Erecs+Etoia;
    iak2=Einfa+EIaenter-EIaleft-Etoia-Ereca;
    ik2i=Einfs+Etoia;
    if2=Einfs+Einfa;
    
    %third step
    Ts2=S(:,tcnt)+sk2/2;
    Te2=E(:,tcnt)+ek2/2;
    Tis2=Is(:,tcnt)+isk2/2;
    Tia2=Ia(:,tcnt)+iak2/2;
    
    ESenter=dt1*(ones(num_loc,1)*theta).*(Migrant*(Ts2./(pop-Tis2-theta*(1-mu).*Tia2)));
    ESleft=min(dt1*(ones(num_loc,1)*theta).*(Ts2./(pop-Tis2-theta*(1-mu).*Tia2)).*(sum(Migrant)'),dt1*Ts2);
    EEenter=dt1*(ones(num_loc,1)*theta.*mu2).*(Migrant*(Te2./(pop-Tis2-theta*(1-mu).*Tia2)));
    EEleft=min(dt1*(ones(num_loc,1)*theta.*mu2).*(Te2./(pop-Tis2-theta*(1-mu).*Tia2)).*(sum(Migrant)'),dt1*Te2);
    EIaenter=dt1*(ones(num_loc,1)*theta.*mu).*(Migrant*(Tia2./(pop-Tis2-theta*(1-mu).*Tia2)));
    EIaleft=min(dt1*(ones(num_loc,1)*theta.*mu).*(Tia2./(pop-Tis2-theta*(1-mu).*Tia2)).*(sum(Migrant)'),dt1*Tia2);
    
    Eexps=dt1*(betas).*Ts2.*Tis2./pop;
    Eexpa=dt1*(betaa).*Ts2.*Tia2./pop;
    Einfs=dt1*(alpha).*Te2./(ones(num_loc,1)*Z);
    Einfa=dt1*((1-alpha)).*Te2./(ones(num_loc,1)*Z);
    Erecs=dt1*Tis2./(D);
    Etoia=dt1*(alpha2).*Tia2;
    Ereca=dt1*Tia2./(D2);
    
    ESenter=max(ESenter,0);ESleft=max(ESleft,0);
    EEenter=max(EEenter,0);EEleft=max(EEleft,0);
    EIaenter=max(EIaenter,0);EIaleft=max(EIaleft,0);
    Eexps=max(Eexps,0);Eexpa=max(Eexpa,0);
    Einfs=max(Einfs,0);Einfa=max(Einfa,0);
    Erecs=max(Erecs,0);Ereca=max(Ereca,0);
    Etoia=max(Etoia,0);
    
    %%%%%%%%%%stochastic version
    ESenter=poissrnd(ESenter);ESleft=poissrnd(ESleft);
    EEenter=poissrnd(EEenter);EEleft=poissrnd(EEleft);
    EIaenter=poissrnd(EIaenter);EIaleft=poissrnd(EIaleft);
    Eexps=poissrnd(Eexps);
    Eexpa=poissrnd(Eexpa);
    Einfs=poissrnd(Einfs);
    Einfa=poissrnd(Einfa);
    Erecs=poissrnd(Erecs);
    Ereca=poissrnd(Ereca);
    Etoia=poissrnd(Etoia);

    sk3=-Eexps-Eexpa+ESenter-ESleft;
    ek3=Eexps+Eexpa-Einfs-Einfa+EEenter-EEleft;
    isk3=Einfs-Erecs+Etoia;
    iak3=Einfa+EIaenter-EIaleft-Etoia-Ereca;
    ik3i=Einfs+Etoia;
    if3=Einfs+Einfa;
    
    %fourth step
    Ts3=S(:,tcnt)+sk3;
    Te3=E(:,tcnt)+ek3;
    Tis3=Is(:,tcnt)+isk3;
    Tia3=Ia(:,tcnt)+iak3;
    
    ESenter=dt1*(ones(num_loc,1)*theta).*(Migrant*(Ts3./(pop-Tis3-theta*(1-mu).*Tia3)));
    ESleft=min(dt1*(ones(num_loc,1)*theta).*(Ts3./(pop-Tis3-theta*(1-mu).*Tia3)).*(sum(Migrant)'),dt1*Ts3);
    EEenter=dt1*(ones(num_loc,1)*theta.*mu2).*(Migrant*(Te3./(pop-Tis3-theta*(1-mu).*Tia3)));
    EEleft=min(dt1*(ones(num_loc,1)*theta.*mu2).*(Te3./(pop-Tis3-theta*(1-mu).*Tia3)).*(sum(Migrant)'),dt1*Te3);
    EIaenter=dt1*(ones(num_loc,1)*theta.*mu).*(Migrant*(Tia3./(pop-Tis3-theta*(1-mu).*Tia3)));
    EIaleft=min(dt1*(ones(num_loc,1)*theta.*mu).*(Tia3./(pop-Tis3-theta*(1-mu).*Tia3)).*(sum(Migrant)'),dt1*Tia3);
    
    Eexps=dt1*(betas).*Ts3.*Tis3./pop;
    Eexpa=dt1*(betaa).*Ts3.*Tia3./pop;
    Einfs=dt1*(alpha).*Te3./(ones(num_loc,1)*Z);
    Einfa=dt1*((1-alpha)).*Te3./(ones(num_loc,1)*Z);
    Erecs=dt1*Tis3./(D);
    Etoia=dt1*(alpha2).*Tia3;
    Ereca=dt1*Tia3./(D2);
    
    ESenter=max(ESenter,0);ESleft=max(ESleft,0);
    EEenter=max(EEenter,0);EEleft=max(EEleft,0);
    EIaenter=max(EIaenter,0);EIaleft=max(EIaleft,0);
    Eexps=max(Eexps,0);Eexpa=max(Eexpa,0);
    Einfs=max(Einfs,0);Einfa=max(Einfa,0);
    Erecs=max(Erecs,0);Ereca=max(Ereca,0);
    Etoia=max(Etoia,0);
    
    %%%%%%%%%%stochastic version
    ESenter=poissrnd(ESenter);ESleft=poissrnd(ESleft);
    EEenter=poissrnd(EEenter);EEleft=poissrnd(EEleft);
    EIaenter=poissrnd(EIaenter);EIaleft=poissrnd(EIaleft);
    Eexps=poissrnd(Eexps);
    Eexpa=poissrnd(Eexpa);
    Einfs=poissrnd(Einfs);
    Einfa=poissrnd(Einfa);
    Erecs=poissrnd(Erecs);
    Ereca=poissrnd(Ereca);
    Etoia=poissrnd(Etoia);


    sk4=-Eexps-Eexpa+ESenter-ESleft;
    ek4=Eexps+Eexpa-Einfs-Einfa+EEenter-EEleft;
    isk4=Einfs-Erecs+Etoia;
    iak4=Einfa+EIaenter-EIaleft-Etoia-Ereca;
    ik4i=Einfs+Etoia;
    if4=Einfs+Einfa;
    
    %%%%%
    S(:,tcnt+1)=S(:,tcnt)+round(sk1/6+sk2/3+sk3/3+sk4/6);
    E(:,tcnt+1)=E(:,tcnt)+round(ek1/6+ek2/3+ek3/3+ek4/6);
    Is(:,tcnt+1)=Is(:,tcnt)+round(isk1/6+isk2/3+isk3/3+isk4/6);
    Ia(:,tcnt+1)=Ia(:,tcnt)+round(iak1/6+iak2/3+iak3/3+iak4/6);
    Incidence(:,tcnt+1)=round(ik1i/6+ik2i/3+ik3i/3+ik4i/6);
    infection(:,tcnt+1)=round(if1/6+if2/3+if3/3+if4/6);
    obs=Incidence(:,tcnt+1);
    
end
%%update pop
pop=pop-sum(Migrant,1)'.*theta+sum(Migrant,2).*theta;
minfrac=0.6;
maxfrac=1.3;
pop(pop<minfrac*pop0)=pop0(pop<minfrac*pop0)*minfrac;
pop(pop>maxfrac*pop0)=pop0(pop>maxfrac*pop0)*maxfrac;
pop2=S(:,2)+E(:,2)+infection(:,2);
k=(pop2>pop);
S(k,tcnt+1)=S(k,tcnt+1).*pop(k)./pop2(k);

%%%update x
x(Sidx)=round(S(:,tcnt+1));
x(Eidx)=round(E(:,tcnt+1));
x(Isidx)=round(Is(:,tcnt+1));
x(Iaidx)=round(Ia(:,tcnt+1));
x(obsidx)=round(sum(obs,2));
x(infidx)=round(infection(:,tcnt+1));
%%%update pop


All=squeeze(E(:,2)+sum(infect,2));
k=(All>pop0);
M(k,:,ts+1:end)=0;
M(:,k,ts+1:end)=0;
end
