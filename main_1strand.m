% Please check line 56 and file gamma_func before running.

clear all
tic

Data=csvread('minimum_data.csv',1 );

Tij=Data(:,2:24);
n=size(Tij,1);
tn=size(Tij,2);
Yij=Data(:,25:47);
Si=Data(:,52);
Sex=Data(:,53);
Age=Data(:,50);
Race1=Data(:,54);
Race2=Data(:,55);
Trt=Data(:,56);
Z1=zeros(n,1);
Z2=zeros(n,1);
Z3=zeros(n,1);
C=zeros(n,1)+10;
for i=1:n
 if  ( Si(i)==0)  Si(i)=3;   C(i)=10; end %
  if (Trt(i)==1) Z1(i)=1; else Z1(i)=0;end %
  if (Trt(i)==6 || Trt(i)==7) C(i)=Si(i); else C(i)=10;end %
  if (Trt(i)==2) Z2(i)=1; else Z2(i)=0;end %
  if (Trt(i)==3) Z3(i)=1; else Z3(i)=0;end %
end  
 
trt1ij=repmat(Z1,1,tn);
trt2ij=repmat(Z2,1,tn);
trt3ij=repmat(Z3,1,tn);
Sexij=repmat(Sex,1,tn);
Race1ij=repmat(Race1,1,tn);
Race2ij=repmat(Race2,1,tn);
Ageij=repmat(Age,1,tn);

Y=sqrt(Yij);
T=Tij;
U=Si;
TU=Tij-repmat(U,1,tn);
n=size(Y,1);
tn=size(Y,2);
X=[ones(n,tn) ]; 
Z=[Sexij Ageij  Race1ij  Race2ij   ]; 
W=[trt1ij.*(TU>0) trt2ij.*(TU>0) trt3ij.*(TU>0)]; 
V=[Z W];
dN=(Y>0).*(T<=repmat(C,1,size(T,2)));
T0=(0.05:0.05:2.45)'; % time points to be estimated 
U0=(0.1:0.05:2.15)';   
N=length(T0);
dt=T0(2:N)'-T0(1:N-1)';
m=sum(dN,2);
h=0.47;%bandwidth
link_id=1; %Set the indication of link function
theta_n=9; %number treatment changing parameters 
alpha_n=size(X,2)/tn; %number of covariates X
beta_n=size(Z,2)/tn;  %number fixed coefficient parameters
gamma_n=size(W,2)/tn; %number 2nd scale varying coefficient parameters
[t1,t2]=make_lim(max(max(T)),h/2);
dN_int=dN.*double((T>t1&T<=t2));
VV=[V,repmat(V(:,tn*beta_n+1:tn*(beta_n+gamma_n)),1,theta_n-gamma_n)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%bandwidth selection%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TT=Tij(:,1);
% TT1=Tij(:,2:23);
% TT1(TT1==0)=nan;
% TTT=[TT,TT1];
% TTTT=TTT;
% VART=nanvar(TTTT,0,2);
% VART(VART==0)=nan;
% MEANT=nanmean(TTTT,2);
% VART1=nanvar(MEANT)+nanmean(VART);
% K=3;%numbers of folds.
% hi=0.35:0.01:0.5;
% [CV,CV1,CV2,CV3,cv_hist]=PE(X,VV,Y,T,U,TU,dN,T0,hi,m,1,link_id,alpha_n,beta_n,theta_n,gamma_n,K,dN_int);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%Main estimation%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=0.41;
alpha_hat0=ones(n,2*alpha_n*tn); %initial values 
beta0=zeros(beta_n+theta_n,1);%initial values 
m2=ones(n,1);
TU=TU.*(TU>0);
[alpha_hat,beta_hat,sig_beta,sig_alpha,alpha_hatT,cov_beta,p_value1,p_value2]=calculate_nowei_test(X,VV,Y,T,U,TU,dN,T0,tn,m,h,alpha_n,n,m2,link_id,beta_n,theta_n,gamma_n,alpha_hat0,beta0,N,U0,dN_int,1000);
[beta_hat_wei,sig_beta_wei,cov_beta_wei]=calculate_wei(X,VV,Y,T,U,TU,dN,T0,tn,m,h,alpha_n,n,ones(n,1),link_id,beta_n,theta_n,gamma_n,alpha_hatT,beta_hat,N,dN_int); 

b_result=[beta_hat,sig_beta, beta_hat-1.96*sig_beta, beta_hat+1.96*sig_beta, 2*(1-normcdf(abs(beta_hat./sig_beta)))]
b_result_wei=[beta_hat_wei,sig_beta_wei, beta_hat_wei-1.96*sig_beta_wei, beta_hat_wei+1.96*sig_beta_wei,2*(1-normcdf(abs(beta_hat_wei./sig_beta_wei)))]

a_result=[alpha_hat,alpha_hat-1.96*sig_alpha,alpha_hat+1.96*sig_alpha];
[resi ,predY]=residual(n,tn,alpha_n,beta_n,theta_n,gamma_n,alpha_hatT,beta_hat ,X,T,U,TU,VV,Y,link_id,dN);
[resi_wei,predY]=residual(n,tn,alpha_n,beta_n,theta_n,gamma_n,alpha_hatT,beta_hat,X,T,U,TU,VV,Y,link_id,dN);
 
N_U=length(U0);

%%%%%%%%%quadratic%%%%%%%%%%%%%%%%%%%%%%%%%
sig_c=U0;c_hat=U0;
for i=1:N_U
 j=U0(i);
 var_c=[1,j,j^2]*cov_beta([beta_n+1,beta_n+theta_n/3+1,beta_n+theta_n/3*2+1], [beta_n+1,beta_n+theta_n/3+1,beta_n+theta_n/3*2+1])*[1;j;j^2];
 sig_c(i)=sqrt(var_c);
 c_hat(i)=beta_hat(beta_n+1)+beta_hat(beta_n+theta_n/3+1)*j+beta_hat(beta_n+theta_n/3*2+1)*j^2;
end
c_result1=[c_hat,c_hat-1.96*sig_c,c_hat+1.96*sig_c];
sig_c=U0;c_hat=U0;
for i=1:N_U
 j=U0(i);
 var_c=[1,j,j^2]*cov_beta([beta_n+2,beta_n+theta_n/3+2,beta_n+theta_n/3*2+2], [beta_n+2,beta_n+theta_n/3+2,beta_n+theta_n/3*2+2])*[1;j;j^2];
 sig_c(i)=sqrt(var_c);
 c_hat(i)=beta_hat(beta_n+2)+beta_hat(beta_n+theta_n/3+2)*j+beta_hat(beta_n+theta_n/3*2+2)*j^2;
end
c_result2=[c_hat,c_hat-1.96*sig_c,c_hat+1.96*sig_c];
sig_c=U0;c_hat=U0;
for i=1:N_U
 j=U0(i);
 var_c=[1,j,j^2]*cov_beta([beta_n+3,beta_n+theta_n/3+3,beta_n+theta_n/3*2+3], [beta_n+3,beta_n+theta_n/3+3,beta_n+theta_n/3*2+3])*[1;j;j^2];
 sig_c(i)=sqrt(var_c);
 c_hat(i)=beta_hat(beta_n+3)+beta_hat(beta_n+theta_n/3+3)*j+beta_hat(beta_n+theta_n/3*2+3)*j^2;
end
c_result3=[c_hat,c_hat-1.96*sig_c,c_hat+1.96*sig_c];


%%%%%%%%%quadratic%%%%%%%%%%%%%%%%%%%%%%%%%
sig_c_wei=U0;c_hat_wei=U0;
for i=1:N_U
 j=U0(i);
 var_c_wei=[1,j,j^2]*cov_beta_wei([beta_n+1,beta_n+theta_n/3+1,beta_n+theta_n/3*2+1], [beta_n+1,beta_n+theta_n/3+1,beta_n+theta_n/3*2+1])*[1;j;j^2];
 sig_c(i)=sqrt(var_c_wei);
 c_hat(i)=beta_hat_wei(beta_n+1)+beta_hat_wei(beta_n+theta_n/3+1)*j+beta_hat_wei(beta_n+theta_n/3*2+1)*j^2;
end
c_result1_wei=[c_hat,c_hat-1.96*sig_c,c_hat+1.96*sig_c];
sig_c_wei=U0;c_hat_wei=U0;
for i=1:N_U
 j=U0(i);
 var_c_wei=[1,j,j^2]*cov_beta_wei([beta_n+2,beta_n+theta_n/3+2,beta_n+theta_n/3*2+2], [beta_n+2,beta_n+theta_n/3+2,beta_n+theta_n/3*2+2])*[1;j;j^2];
 sig_c(i)=sqrt(var_c_wei);
 c_hat(i)=beta_hat_wei(beta_n+2)+beta_hat_wei(beta_n+theta_n/3+2)*j+beta_hat_wei(beta_n+theta_n/3*2+2)*j^2;
end
c_result2_wei=[c_hat,c_hat-1.96*sig_c,c_hat+1.96*sig_c];
sig_c_wei=U0;c_hat_wei=U0;
for i=1:N_U
 j=U0(i);
 var_c_wei=[1,j,j^2]*cov_beta_wei([beta_n+3,beta_n+theta_n/3+3,beta_n+theta_n/3*2+3], [beta_n+3,beta_n+theta_n/3+3,beta_n+theta_n/3*2+3])*[1;j;j^2];
 sig_c(i)=sqrt(var_c_wei);
 c_hat(i)=beta_hat_wei(beta_n+3)+beta_hat_wei(beta_n+theta_n/3+3)*j+beta_hat_wei(beta_n+theta_n/3*2+3)*j^2;
end
c_result3_wei=[c_hat,c_hat-1.96*sig_c,c_hat+1.96*sig_c];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc  

% plot(T0,a_result((1:N),1),T0,a_result((1:N),2),T0,a_result((1:N),3))
% plot(T0,a_result((1:N),4),T0,a_result((1:N),5),T0,a_result((1:N),6))
% plot(Y(resi~=0),resi(resi~=0),'.')

