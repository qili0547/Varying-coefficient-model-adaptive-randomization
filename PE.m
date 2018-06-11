function [CV,CV1,CV2,CV3,cv_hist]=PE(X,VV,Y,T,U,TU,dN,T0,hi,m,rnd, link_id,alpha_n,beta_n,theta_n,gamma_n,K,dN_int)

%This function is used to calculate the CV values.
%rnd is the indivator of whether use a random arrangement for the folds.
%K is the number of folds.

tn=size(T,2);
n=size(T,1);
N=length(T0);
cv=zeros(K,1);CV=zeros(size(hi));
beta_hist=zeros(length(hi),K*beta_n);
L0=(n/K);
if rnd==0 [a,id]=sort(1:n);
else [a,id]=sort(rand(n,1)); 
end
cv_hist=zeros(length(hi),K);
for i=1:length(hi)
    h1=hi(i);
    for k=1:K
          id1=id;
          id2=fix((k-1)*L0)+1:fix(k*L0);
          X2=X(id1(id2),:);  
          V2=VV(id1(id2),:); 
          Y2=Y(id1(id2),:); 
          T2=T(id1(id2),:);
		  U2=U(id1(id2),:);
		  TU2=TU(id1(id2),:);
          dN2=dN(id1(id2),:); 
		  dN_int2=dN_int(id1(id2),:); 
          m2=m(id1(id2),:); 
		  n2=size(X2,1);
          %%%%%%%%
           id1(id2)=[];
          X1=X(id1,:);  %Leave k-group out%
          V1=VV(id1,:); 
          Y1=Y(id1,:); 
          T1=T(id1,:); 
		  U1=U(id1,:);
		  TU1=TU(id1,:);
          dN1=dN(id1,:); 
		  dN_int1=dN_int(id1,:); 
          m1=m(id1,:); 
          n1=size(X1,1);
		  
     alpha_hat0=ones(n1,2*alpha_n*tn); %initial values%
     beta0=zeros(beta_n+theta_n,1);%initial values%
    [alpha_hat,beta_hat]=calculate_nowei(X1,V1,Y1,T1,U1,TU1,dN1,T2,tn,m1,h1,alpha_n,n1,ones(n1,1),link_id,beta_n,theta_n,gamma_n,alpha_hat0,beta0,N,dN_int1);
    %beta_hist(i,k:K:end)=beta_hat;
    
	temp=0; 
    for id_gamma=1:alpha_n
       temp= temp+alpha_hat(:,tn*(id_gamma-1)+(1:tn)).*X2(:,tn*(id_gamma-1)+(1:tn));
    end
	zeta001=zeros(n2,tn*beta_n);
	for id_beta=1:beta_n
		   zeta001(1:n2,tn*(id_beta-1)+(1:tn))=repmat(beta0(id_beta),n2,tn);
	end
	
	[zeta002,dzeta002]=gamma_func(beta_hat(beta_n+1:beta_n+theta_n),T2,U2,TU2,theta_n,gamma_n);
	if beta_n>0  zeta012=[zeta001,zeta002]; else zeta012=[zeta002]; end
	temp02=0;
	for id_beta=1:beta_n+theta_n;
	  temp02= temp02+zeta012(:,tn*(id_beta-1)+(1:tn)).* V2(:,tn*(id_beta-1)+(1:tn));
	end
    [mu_hat,mu_dot_hat]=psy(temp+temp02,link_id);
    cv(k)=sum(sum((Y2-mu_hat).^2.*dN2));
    cv1(k)=sum(sum(mu_dot_hat.*(Y2-mu_hat).^2.*dN2));
    cv2(k)=sum(sum((Y2-mu_hat).*dN2).^2);
    cv3(k)=sum(sum(mu_dot_hat.*(Y2-mu_hat).*dN2).^2);
    end
     CV(i)=sum(cv);
	 CV1(i)=sum(cv1);
	 CV2(i)=sum(cv2);
	 CV3(i)=sum(cv3);
     cv_hist(i,:)=cv;
end
 