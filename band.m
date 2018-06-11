clear all
tic

Data=csvread('wide1.csv',1,0);

T=Data(:,8:19);
Y=Data(:,20:31);
Xcd=Data(:,32:43);
Zad5=log(Data(:,44:55));
Zad5(~isfinite(Zad5))=0;
Xtrt=Data(:,56:67);
Zaustnam=Data(:,68:79);
Vii=Data(:,80:91);
Cii=Data(:,104:115);
Sii=Data(:,116:127);
ageii=Data(:,128:139);
lastvacii=Data(:,140:151);
ppii=Data(:,152:163);

Y=Y;
T=T;
n=size(Y,1);
tn=size(Y,2);
X=[ones(n,tn) Xcd];
Z=[Zad5  Zaustnam ppii];
W=[Xtrt];
V=[Z W];
U=-ageii;
dN=(ageii>0)*1;
T0=(0.05:0.05:2.5)';  % time points to be estimated
N=length(T0);
dt=T0(2:N)'-T0(1:N-1)';
m=sum(dN,2);
h=.5;

alpha_n=size(X,2)/tn; %number of covariates X
beta_n=size(Z,2)/tn;  %number fixed coefficient parameters
gamma_n=size(W,2)/tn; %number 2nd scale varying coefficient parameters
theta_n=2; %number treatment changing parameters

link_id=1; %Set the indication of link function(see more in the function psy.m)
VV=[V,repmat(V(:,tn*beta_n+1:tn*(beta_n+gamma_n)),1,1)];
[t1,t2]=make_lim(max(max(T)),h);
dN_int=dN;

K=10;%K is the set of numbers of folds.
hi=0.2:0.01:0.8;
[CV,CV1,CV2,CV3,cv_hist]=PE(X,VV,Y,T,U,dN,T0,hi,m,1,link_id,alpha_n,beta_n,theta_n,gamma_n,K,dN_int);

save example1_band.mat CV CV1 CV2 CV3 ;

toc  
