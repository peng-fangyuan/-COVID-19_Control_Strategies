function [x,pop]=SEIR(x,M,innerflow,pop,ts,pop0)
%the metapopulation SEIR model
load hu
dt=1;
tmstep=1;
%integrate forward for one day
num_loc=size(pop,1);
[~,num_ens]=size(x);
%S,E,Is,Ia,obs,...,beta,mu,Z,alpha,D
Sidx=(1:5:5*num_loc)';
Eidx=(2:5:5*num_loc)';
Isidx=(3:5:5*num_loc)';
Iaidx=(4:5:5*num_loc)';
obsidx=(5:5:5*num_loc)';
betaidx=5*num_loc+2;
muidx=5*num_loc+1;
muiidx=5*num_loc+3;
Zidx=5*num_loc+4;
alphaidx=5*num_loc+5;
Didx=5*num_loc+6;
D2idx=5*num_loc+7;

alphaidx2=5*num_loc+8; % testing rate for wuhan
S=zeros(num_loc,num_ens,tmstep+1);
E=zeros(num_loc,num_ens,tmstep+1);
Is=zeros(num_loc,num_ens,tmstep+1);
Ia=zeros(num_loc,num_ens,tmstep+1);
Incidence=zeros(num_loc,num_ens,tmstep+1);
obs=zeros(num_loc,num_ens);
%initialize S,E,Is,Ia and parameters
S(:,:,1)=x(Sidx,:);
E(:,:,1)=x(Eidx,:);
Is(:,:,1)=x(Isidx,:);
Ia(:,:,1)=x(Iaidx,:);


theta=ones(num_loc,1)*6.5*ones(1,num_ens);
mui=x(muiidx,:);
Z=x(Zidx,:);
alpha=ones(num_loc,1)*x(alphaidx,:);
D=ones(num_loc,1)*x(Didx,:);
D2=ones(num_loc,1)*x(D2idx,:);
mu=x(muidx,:);
report=ones(num_loc,1)*x(alphaidx,:); %incidence in previous days being tested possibility
alpha(292,:)=alpha(292,:).*x(alphaidx2,:); %testing rate for wuhan
report(292,:)=report(292,:).*x(alphaidx2,:);

if ts==34 % model the scenerio that in 2.13, cases are confirmed based on clinical diganosis cities in Huber and the number sharply increase
    for i=1:length(hu)
        loc=hu(i);
        report(loc,:)=1;
    end
end


%start integration
tcnt=0;

betaa=innerflow(:,ts)*x(betaidx,:);
betas=betaa.*mu;
Migrant=M(:,:,ts);
for t=ts+dt:dt:ts+tmstep
    tcnt=tcnt+1;
    dt1=dt;
    %first step
    ESenter=dt1*(theta).*(Migrant*(S(:,:,tcnt)./(pop-Is(:,:,tcnt))));
    ESleft=min(dt1*(theta).*(S(:,:,tcnt)./(pop-Is(:,:,tcnt))).*(sum(Migrant)'*ones(1,num_ens)),dt1*S(:,:,tcnt));
    EEenter=dt1*(theta.*mui).*(Migrant*(E(:,:,tcnt)./(pop-Is(:,:,tcnt))));
    EEleft=min(dt1*(theta.*mui).*(E(:,:,tcnt)./(pop-Is(:,:,tcnt))).*(sum(Migrant)'*ones(1,num_ens)),dt1*E(:,:,tcnt));
    EIaenter=dt1*(theta.*mui).*(Migrant*(Ia(:,:,tcnt)./(pop-Is(:,:,tcnt))));
    EIaleft=min(dt1*(theta.*mui).*(Ia(:,:,tcnt)./(pop-Is(:,:,tcnt))).*(sum(Migrant)'*ones(1,num_ens)),dt1*Ia(:,:,tcnt));
    
    Eexps=dt1*(betas).*S(:,:,tcnt).*Is(:,:,tcnt)./pop;
    Eexpa=dt1*(betaa).*S(:,:,tcnt).*Ia(:,:,tcnt)./pop;
    Einfs=dt1*(alpha).*E(:,:,tcnt)./(ones(num_loc,1)*Z);
    Einfa=dt1*((1-alpha)).*E(:,:,tcnt)./(ones(num_loc,1)*Z);
    Erecs=dt1*Is(:,:,tcnt)./(D);
    Etoia=dt1*(report).*Ia(:,:,tcnt);
    Ereca=dt1*Ia(:,:,tcnt)./(D2);
    
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
    iak1=Einfa+EIaenter-EIaleft-Ereca-Etoia;
    ik1i=Einfs+Etoia;
    %second step
    Ts1=S(:,:,tcnt)+sk1/2;
    Te1=E(:,:,tcnt)+ek1/2;
    Tis1=Is(:,:,tcnt)+isk1/2;
    Tia1=Ia(:,:,tcnt)+iak1/2;
    
    ESenter=dt1*(theta).*(Migrant*(Ts1./(pop-Tis1)));
    ESleft=min(dt1*(theta).*(Ts1./(pop-Tis1)).*(sum(Migrant)'*ones(1,num_ens)),dt1*Ts1);
    EEenter=dt1*(theta.*mui).*(Migrant*(Te1./(pop-Tis1)));
    EEleft=min(dt1*(theta.*mui).*(Te1./(pop-Tis1)).*(sum(Migrant)'*ones(1,num_ens)),dt1*Te1);
    EIaenter=dt1*(theta.*mui).*(Migrant*(Tia1./(pop-Tis1)));
    EIaleft=min(dt1*(theta.*mui).*(Tia1./(pop-Tis1)).*(sum(Migrant)'*ones(1,num_ens)),dt1*Tia1);
    
    Eexps=dt1*(betas).*Ts1.*Tis1./pop;
    Eexpa=dt1*(betaa).*Ts1.*Tia1./pop;
    Einfs=dt1*(alpha).*Te1./(ones(num_loc,1)*Z);
    Einfa=dt1*((1-alpha)).*Te1./(ones(num_loc,1)*Z);
    Erecs=dt1*Tis1./(D);
    Etoia=dt1*(report).*Tia1;
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
    iak2=Einfa+EIaenter-EIaleft-Ereca-Etoia;
    ik2i=Einfs+Etoia;
    
    %third step
    Ts2=S(:,:,tcnt)+sk2/2;
    Te2=E(:,:,tcnt)+ek2/2;
    Tis2=Is(:,:,tcnt)+isk2/2;
    Tia2=Ia(:,:,tcnt)+iak2/2;
    
    ESenter=dt1*(theta).*(Migrant*(Ts2./(pop-Tis2)));
    ESleft=min(dt1*(theta).*(Ts2./(pop-Tis2)).*(sum(Migrant)'*ones(1,num_ens)),dt1*Ts2);
    EEenter=dt1*(theta.*mui).*(Migrant*(Te2./(pop-Tis2)));
    EEleft=min(dt1*(theta.*mui).*(Te2./(pop-Tis2)).*(sum(Migrant)'*ones(1,num_ens)),dt1*Te2);
    EIaenter=dt1*(theta.*mui).*(Migrant*(Tia2./(pop-Tis2)));
    EIaleft=min(dt1*(theta.*mui).*(Tia2./(pop-Tis2)).*(sum(Migrant)'*ones(1,num_ens)),dt1*Tia2);
    
    Eexps=dt1*(betas).*Ts2.*Tis2./pop;
    Eexpa=dt1*(betaa).*Ts2.*Tia2./pop;
    Einfs=dt1*(alpha).*Te2./(ones(num_loc,1)*Z);
    Einfa=dt1*((1-alpha)).*Te2./(ones(num_loc,1)*Z);
    Erecs=dt1*Tis2./(D);
    Etoia=dt1*(report).*Tia2;
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
    iak3=Einfa+EIaenter-EIaleft-Ereca-Etoia;
    ik3i=Einfs+Etoia;
    
    %fourth step
    Ts3=S(:,:,tcnt)+sk3;
    Te3=E(:,:,tcnt)+ek3;
    Tis3=Is(:,:,tcnt)+isk3;
    Tia3=Ia(:,:,tcnt)+iak3;
    
    ESenter=dt1*(theta).*(Migrant*(Ts3./(pop-Tis3)));
    ESleft=min(dt1*(theta).*(Ts3./(pop-Tis3)).*(sum(Migrant)'*ones(1,num_ens)),dt1*Ts3);
    EEenter=dt1*(theta.*mui).*(Migrant*(Te3./(pop-Tis3)));
    EEleft=min(dt1*(theta.*mui).*(Te3./(pop-Tis3)).*(sum(Migrant)'*ones(1,num_ens)),dt1*Te3);
    EIaenter=dt1*(theta.*mui).*(Migrant*(Tia3./(pop-Tis3)));
    EIaleft=min(dt1*(theta.*mui).*(Tia3./(pop-Tis3)).*(sum(Migrant)'*ones(1,num_ens)),dt1*Tia3);
    
    Eexps=dt1*(betas).*Ts3.*Tis3./pop;
    Eexpa=dt1*(betaa).*Ts3.*Tia3./pop;
    Einfs=dt1*(alpha).*Te3./(ones(num_loc,1)*Z);
    Einfa=dt1*((1-alpha)).*Te3./(ones(num_loc,1)*Z);
    Erecs=dt1*Tis3./(D);
    Etoia=dt1*(report).*Tia3;
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
    iak4=Einfa+EIaenter-EIaleft-Ereca-Etoia;
    ik4i=Einfs+Etoia;
    
    %%%%%
    S(:,:,tcnt+1)=S(:,:,tcnt)+round(sk1/6+sk2/3+sk3/3+sk4/6);
    E(:,:,tcnt+1)=E(:,:,tcnt)+round(ek1/6+ek2/3+ek3/3+ek4/6);
    Is(:,:,tcnt+1)=Is(:,:,tcnt)+round(isk1/6+isk2/3+isk3/3+isk4/6);
    Ia(:,:,tcnt+1)=Ia(:,:,tcnt)+round(iak1/6+iak2/3+iak3/3+iak4/6);
    Incidence(:,:,tcnt+1)=round(ik1i/6+ik2i/3+ik3i/3+ik4i/6);
    obs=Incidence(:,:,tcnt+1);
end
%%%update x
x(Sidx,:)=S(:,:,tcnt+1);
x(Eidx,:)=E(:,:,tcnt+1);
x(Isidx,:)=Is(:,:,tcnt+1);
x(Iaidx,:)=Ia(:,:,tcnt+1);
x(obsidx,:)=obs;
%%%update pop
pop=pop-sum(M(:,:,ts),1)'.*theta+sum(M(:,:,ts),2).*theta;
minfrac=0.6;
pop(pop<minfrac*pop0)=pop0(pop<minfrac*pop0)*minfrac;