function [alpha_hat,beta_hat,sig_beta,sig_alpha,alpha_hatT,cov_beta,p_ks1,p_ks2]=Calculate(X,VV,Y,T,U,TU,dN,T0,tn,m,h,alpha_n,n,m2,link_id,beta_n,theta_n,gamma_n,alpha_hat0,beta0,N,U0,dN_int,D)

N0=sum(m);
alpha_hat1=alpha_hat0;
epsilon=0.001;
id=1;
while(id<20)
 
   for(i=1:n)
   for(j=1:m(i))
      alpha0=alpha_hat0(i,j:tn:(2*alpha_n*tn))';
	  t=T(i,j);
	  alpha1=alpha0;
	  iter=1; Error=1;
      if sum(sum(T>t-h&T<t+h&T>0))>=5
      while( Error>epsilon & iter < 50)
         temp00=0; temp01=0;temp02=0;
         for id_alpha=1:alpha_n
            temp00= temp00+alpha0(id_alpha)*X(:,tn*(id_alpha-1)+(1:tn));
            temp01= temp01+alpha0(alpha_n+id_alpha)*X(:,tn*(id_alpha-1)+(1:tn)).*(T-t);
         end
		 for id_beta=1:beta_n
		   zeta01(1:n,tn*(id_beta-1)+(1:tn))=repmat(beta0(id_beta),n,tn);
		 end
		 [zeta02,dzeta02]=gamma_func(beta0(beta_n+1:beta_n+theta_n),T,U,TU,theta_n,gamma_n);
		 if beta_n>0  zeta012=[zeta01,zeta02]; else zeta012=[zeta02];end
		 for id_beta=1:beta_n+theta_n
		    temp02=temp02+zeta012(:,tn*(id_beta-1)+(1:tn)).*VV(:,tn*(id_beta-1)+(1:tn));     
         end		 
		 [mu_a,mu_dot]=psy(temp00+temp01+temp02,link_id);
         [Ua,Jac]=Ualpha((Y-mu_a).*dN,[X,X.*(repmat(T,1,alpha_n)-t)],T,t,h,alpha_n,tn,m,mu_dot.*dN);
         if (rcond(Jac)<0.00000001)    break; end
         try  alpha1=alpha0-pinv(Jac)*Ua;
              Error=sum(abs(alpha1-alpha0));iter=iter+1;
              alpha0=alpha1;
              catch break;
              end
      end
      end
      alpha_hat1(i,j:tn:(2*alpha_n*tn))=alpha1';
    end
    end

    alpha=alpha_hat1;
    iter=1;Error=1;beta00=beta0;beta11=beta00;
    temp=0; 
    for id_alpha=1:alpha_n
       temp= temp+alpha(:,tn*(id_alpha-1)+(1:tn)).*X(:,tn*(id_alpha-1)+(1:tn));
    end
    while( Error>epsilon & iter < 20)
	   	
		for id_beta=1:beta_n
		   zeta001(1:n,tn*(id_beta-1)+(1:tn))=repmat(beta00(id_beta),n,tn);
		end
	   [zeta002,dzeta002]=gamma_func(beta00(beta_n+1:beta_n+theta_n),T,U,TU,theta_n,gamma_n);
	   
	   if beta_n>0  zeta00=[zeta001,zeta002]; else zeta00=[zeta002];end
	   dzeta001=repmat(1,n,beta_n*tn);
	   dzeta00=[dzeta001,dzeta002];
	   [Pgx,IExx,Exz]=Pgamma_b_X(alpha,X,VV,Y,T,dN,tn,m,h,zeta00,dzeta00,T0,0,link_id);
	   %Pgx is dalpha/dgamma*X(t).
	   temp02=0;
	   for id_beta=1:beta_n+theta_n
	    temp02= temp02+zeta00(:,tn*(id_beta-1)+(1:tn)).* VV(:,tn*(id_beta-1)+(1:tn));
	   end
       [mu,mu_dot_beta]=psy(temp+temp02,link_id);
	   
       Ubeta=zeros(beta_n+theta_n,1);     Jacb=zeros(beta_n+theta_n,beta_n+theta_n);
       for id_gamma=1:(beta_n+theta_n)
           Ubeta(id_gamma)=sum(sum((Y-mu).*(Pgx(:,tn*(id_gamma-1)+(1:tn))+dzeta00(:,tn*(id_gamma-1)+(1:tn)).*VV(:,tn*(id_gamma-1)+(1:tn))).*dN_int));
	      for id_beta=1:(beta_n+theta_n)
               temp02=(Pgx(:,tn*(id_gamma-1)+(1:tn))+dzeta00(:,tn*(id_gamma-1)+(1:tn)).*VV(:,tn*(id_gamma-1)+(1:tn))).*(Pgx(:,tn*(id_beta-1)+(1:tn))+dzeta00(:,tn*(id_beta-1)+(1:tn)).*VV(:,tn*(id_beta-1)+(1:tn)));
               Jacb(id_gamma,id_beta)=-sum(sum(mu_dot_beta.*temp02.*dN_int));
           end
       end
       if (rcond(Jacb)<0.00000001)    break;   end
       try  beta11=beta00-pinv(Jacb)*Ubeta;
            Error=sum(abs(beta11-beta00));iter=iter+1;
            beta00=beta11;
       catch break;
       end
     end
	beta1=beta11;
	
    MISE=sum(sum(abs(alpha_hat1-alpha_hat0).*repmat(dN,1,2*alpha_n)))/N0+norm(beta1-beta0);
    if (MISE>0.005)
        beta0=beta1;
        alpha_hat0=alpha_hat1;
        id=id+1;        
     else break;
    end            
end
beta_hat=beta1;

if nargout>2
%%%%% Estimating the variances of beta
    temp=0; temp02=0; 
    for id_alpha=1:alpha_n
       temp= temp+alpha(:,tn*(id_alpha-1)+(1:tn)).*X(:,tn*(id_alpha-1)+(1:tn));
    end
    for id_beta=1:beta_n
	   zeta001(1:n,tn*(id_beta-1)+(1:tn))=repmat(beta_hat(id_beta),n,tn);
	end
	dzeta001=repmat(1,n,tn*beta_n);
	[zeta002,dzeta002]=gamma_func(beta_hat(beta_n+1:beta_n+theta_n),T,U,TU,theta_n,gamma_n);
	if beta_n>0  zeta00=[zeta001,zeta002]; else zeta00=[zeta002];end
	dzeta00=[dzeta001,dzeta002];
	temp02=0;
	for id_beta=1:beta_n+theta_n
	  temp02= temp02+zeta00(:,tn*(id_beta-1)+(1:tn)).* VV(:,tn*(id_beta-1)+(1:tn));
	end
	[Pgx,IExx,Exz]=Pgamma_b_X(alpha,X,VV,Y,T,dN,tn,m,h,zeta00,dzeta00,T0,0,link_id);
    [mu_hat,mu_dot_beta]=psy(temp+temp02,link_id);
    A=zeros(beta_n+theta_n,beta_n+theta_n);Sigma=A;
    for id_gamma=1:(beta_n+theta_n)
    for id_beta=1:(beta_n+theta_n)
        temp02=(Pgx(:,tn*(id_gamma-1)+(1:tn))+dzeta00(:,tn*(id_gamma-1)+(1:tn)).*VV(:,tn*(id_gamma-1)+(1:tn))).*(Pgx(:,tn*(id_beta-1)+(1:tn))+dzeta00(:,tn*(id_beta-1)+(1:tn)).*VV(:,tn*(id_beta-1)+(1:tn)));
        A(id_gamma,id_beta)=sum(sum(mu_dot_beta.*temp02.*dN_int));
        Sigma(id_gamma,id_beta)=sum(((Y-mu_hat).*(Pgx(:,tn*(id_gamma-1)+(1:tn))+dzeta00(:,tn*(id_gamma-1)+(1:tn)).*VV(:,tn*(id_gamma-1)+(1:tn))).*dN_int),2)'*...
            sum(((Y-mu_hat).*(Pgx(:,tn*(id_beta-1)+(1:tn))+dzeta00(:,tn*(id_beta-1)+(1:tn)).*VV(:,tn*(id_beta-1)+(1:tn))).*dN_int),2);
    end
    end
   sig_beta=diag(sqrt(pinv(A)*Sigma*pinv(A)));
   cov_beta= pinv(A)*Sigma*pinv(A) ;
   
   %%%for the variance of alpha(t)
   sig_alpha=zeros(size(T0,1),alpha_n);
   InvA=pinv(A/n);
   [Pgx1,IExx,Exz]=Pgamma_b_X(alpha_hat1,X,VV,Y,T,dN,tn,m,h,zeta00,dzeta00,T0,1,link_id);
   temp0=zeros(n,beta_n+theta_n);
   for id_beta=1:beta_n+theta_n
       temp0(:,id_beta)= sum(((Y-mu_hat).*(Pgx(:,tn*(id_beta-1)+(1:tn))+dzeta00(:,tn*(id_beta-1)+(1:tn)).*VV(:,tn*(id_beta-1)+(1:tn))).*dN_int),2);
   end
   
end
%Estimate alpha(t) on time points of T0
alpha_hat00=Ker(alpha_hat1,T,T0,alpha_n,tn,h,m,m2); 
tn1=size(T0,2);
alpha_hat11=alpha_hat00;
 for(i=1:size(T0,1))
 for(j=1:m2(i))
      alpha0=alpha_hat00(i,j:tn1:(2*alpha_n*tn1))';
	  t=T0(i,j);
	  alpha1=alpha0;
	  iter=1; Error=1;
      if sum(sum(T>t-h&T<t+h&T>0))>=5
      while( Error>epsilon & iter < 50)
           temp00=0; temp01=0;
           for id_gamma=1:alpha_n
               temp00= temp00+alpha0(id_gamma)*X(:,tn*(id_gamma-1)+(1:tn));
               temp01= temp01+alpha0(alpha_n+id_gamma)*X(:,tn*(id_gamma-1)+(1:tn)).*(T-t);
           end
           temp02=0;

          for id_beta=1:beta_n
	          zeta_hat1(1:n,tn*(id_beta-1)+(1:tn))=repmat(beta_hat(id_beta),n,tn);
	      end
		   
		   [zeta_hat2,dzeta_hat2]=gamma_func(beta_hat(beta_n+1:beta_n+theta_n),T,U,TU,theta_n,gamma_n);
		   if beta_n>0  zeta_hat12=[zeta_hat1,zeta_hat2]; else zeta_hat12=[zeta_hat2]; end
		   for id_beta=1:beta_n+theta_n;
		       temp02= temp02+zeta_hat12(:,tn*(id_beta-1)+(1:tn)).*VV(:,tn*(id_beta-1)+(1:tn));
	       end
           [mu_a,mu_dot]=psy(temp00+temp01+temp02,link_id);
           [Ua,Jac]=Ualpha((Y-mu_a).*dN,[X,X.*(repmat(T,1,alpha_n)-t)],T,t,h,alpha_n,tn,m,mu_dot.*dN);
           if (rcond(Jac)<0.00000001)    break; end
               try  alpha1=alpha0-pinv(Jac)*Ua;
                 Error=sum(abs(alpha1-alpha0));iter=iter+1;
                 alpha0=alpha1;
               catch break;
               end
      end
      end
      alpha_hat11(i,j:tn1:(2*alpha_n*tn1))=alpha1';
	  if nargout>2
         var_alpha=alpha_var(Y,mu_a,X,T,dN,h,t,alpha_n,tn,m,IExx(i,:),Exz(i,:),InvA,temp0*InvA);
         sig_alpha(i,:)=sqrt(diag(var_alpha))/n;
      end
  end
  end
	alpha_hat=alpha_hat11(:,1:alpha_n*tn1);   
%alpha_hat=min(max(alpha_hat,-10),10);    
alpha_hatT=alpha_hat1;
    
	weight=T-T+1;;
    Rij=zeros(size(U0,1),beta_n+theta_n);
    randn('seed',110);
	G=normrnd(0,1,n,D);
 
	
	A=zeros(beta_n+theta_n,beta_n+theta_n); 
    for id_gamma=1:(beta_n+theta_n)
    for id_beta=1:(beta_n+theta_n)
        temp02=(Pgx(:,tn*(id_gamma-1)+(1:tn))+dzeta00(:,tn*(id_gamma-1)+(1:tn)).*VV(:,tn*(id_gamma-1)+(1:tn))).*(Pgx(:,tn*(id_beta-1)+(1:tn))+dzeta00(:,tn*(id_beta-1)+(1:tn)).*VV(:,tn*(id_beta-1)+(1:tn)));
        A(id_gamma,id_beta)=nansum(nansum(mu_dot_beta.*temp02.*dN_int))/n;
    end
    end
	
	U0=(0.1:0.1:2.25)';
	%U0=(0.1:0.1:0.7)';
	
  for(i=1:size(U0,1))	
    u=U0(i);
	temp=0;
    for id_gamma=1:(beta_n+theta_n)
	   for id_beta=1:(beta_n+theta_n)
	      temp=(Pgx(:,tn*(id_gamma-1)+(1:tn))+dzeta00(:,tn*(id_gamma-1)+(1:tn)).*VV(:,tn*(id_gamma-1)+(1:tn))).*(Pgx(:,tn*(id_beta-1)+(1:tn))+dzeta00(:,tn*(id_beta-1)+(1:tn)).*VV(:,tn*(id_beta-1)+(1:tn)));
          Au(id_gamma,id_beta)=nansum(nansum(mu_dot_beta.*temp.*(TU<=u).*dN_int))/n;
	   end
       Qij=(Y-mu_hat).*(Pgx(:,tn*(id_gamma-1)+(1:tn))+dzeta00(:,tn*(id_gamma-1)+(1:tn)).*VV(:,tn*(id_gamma-1)+(1:tn)));     
	   Rij0(id_gamma,1)=sum(sum(Qij.*(TU<=u).*dN_int))/sqrt(n);	   
	   Dij1=repmat(sum(Qij.*(TU<=u).*dN_int,2),1,D).*G;
	   Dij2=repmat(sum(weight.*Qij.*dN_int,2),1,D).*G;
	   Dij2a(id_gamma,:)=sum(Dij2);
	   Dij1a(id_gamma,:)=sum(Dij1);
    end
	Rij(i,:)=pinv(A)*Rij0;
	Rij_star1(i,:,:)=pinv(A)*Dij1a/sqrt(n);
	Rij_star2(i,:,:)=-pinv(A)*Au*pinv(A)*Dij2a/sqrt(n);
	
  end
   Rij_star=Rij_star1+Rij_star2;

  test1=max(sqrt(sum(Rij.^2,2)));
  test1star=max(sqrt(sum(Rij_star(:,beta_n+1:beta_n+theta_n,:).^2,2)));
  
  test2=sum(sqrt(sum(Rij.^2,2)));
  test2star=sum(sqrt(sum(Rij_star.^2,2)));
  
  p_ks1=mean(test1<test1star);
  p_ks2=mean(test2<test2star);
  
  


   


