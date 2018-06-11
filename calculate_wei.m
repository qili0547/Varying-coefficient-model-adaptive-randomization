function [beta_hat,sig_beta,cov_beta]=Calculate_wei(X,VV,Y,T,U,TU,dN,T0,tn,m,h,alpha_n,n,m2,link_id,beta_n,theta_n,gamma_n,alpha_hat0,beta0,N,dN_int)

N0=sum(m);
alpha_hat1=alpha_hat0;
epsilon=0.001;
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
	   for id_beta=1:beta_n+theta_n;
	    temp02= temp02+zeta00(:,tn*(id_beta-1)+(1:tn)).* VV(:,tn*(id_beta-1)+(1:tn));
	   end
       [mu_beta,mu_dot_beta]=psy(temp+temp02,link_id);	   
	   sigerr=Ker2((Y-mu_beta).^2,T,T,1,tn,h,m,m); 
	   weight=mu_dot_beta./sigerr;
	   weight(isinf(weight)) = 0 ;   
	   if (link_id==2)
	     sigerr=Ker2(mu_dot_beta,T,T,1,tn,h,m,m); 
	     weight=mu_dot_beta./sigerr;
	     weight(isinf(weight)) = 0 ;
	   end
       Ubeta=zeros(beta_n+theta_n,1);     Jacb=zeros(beta_n+theta_n,beta_n+theta_n);
       for id_gamma=1:(beta_n+theta_n)
           Ubeta(id_gamma)=sum(sum(weight.*(Y-mu_beta).*(Pgx(:,tn*(id_gamma-1)+(1:tn))+dzeta00(:,tn*(id_gamma-1)+(1:tn)).*VV(:,tn*(id_gamma-1)+(1:tn))).*dN_int));
	      for id_beta=1:(beta_n+theta_n)
               temp02=(Pgx(:,tn*(id_gamma-1)+(1:tn))+dzeta00(:,tn*(id_gamma-1)+(1:tn)).*VV(:,tn*(id_gamma-1)+(1:tn))).*(Pgx(:,tn*(id_beta-1)+(1:tn))+dzeta00(:,tn*(id_beta-1)+(1:tn)).*VV(:,tn*(id_beta-1)+(1:tn)));
               Jacb(id_gamma,id_beta)=-sum(sum(weight.*mu_dot_beta.*temp02.*dN_int));
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
	iter
beta_hat=beta1;

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
	for id_beta=1:beta_n+theta_n;
	  temp02= temp02+zeta00(:,tn*(id_beta-1)+(1:tn)).* VV(:,tn*(id_beta-1)+(1:tn));
	end
	[Pgx,IExx,Exz]=Pgamma_b_X(alpha,X,VV,Y,T,dN,tn,m,h,zeta00,dzeta00,T0,0,link_id);
    [mu_hat,mu_dot_beta]=psy(temp+temp02,link_id);
      sigerr=Ker2((Y-mu_hat).^2,T,T,1,tn,h,m,m); 
	   weight=mu_dot_beta./sigerr;
	   weight(isinf(weight)) = 0 ;
	if (link_id==2)
	     sigerr=Ker2(mu_dot_beta,T,T,1,tn,h,m,m); 
	     weight=mu_dot_beta./sigerr;
	     weight(isinf(weight)) = 0 ;
	end
	
    A=zeros(beta_n+theta_n,beta_n+theta_n);Sigma=A;
    for id_gamma=1:(beta_n+theta_n)
    for id_beta=1:(beta_n+theta_n)
        temp02=(Pgx(:,tn*(id_gamma-1)+(1:tn))+dzeta00(:,tn*(id_gamma-1)+(1:tn)).*VV(:,tn*(id_gamma-1)+(1:tn))).*(Pgx(:,tn*(id_beta-1)+(1:tn))+dzeta00(:,tn*(id_beta-1)+(1:tn)).*VV(:,tn*(id_beta-1)+(1:tn)));
        A(id_gamma,id_beta)=sum(sum(weight.*mu_dot_beta.*temp02.*dN_int));
        Sigma(id_gamma,id_beta)=sum((weight.*(Y-mu_hat).*(Pgx(:,tn*(id_gamma-1)+(1:tn))+dzeta00(:,tn*(id_gamma-1)+(1:tn)).*VV(:,tn*(id_gamma-1)+(1:tn))).*dN_int),2)'*...
            sum((weight.*(Y-mu_hat).*(Pgx(:,tn*(id_beta-1)+(1:tn))+dzeta00(:,tn*(id_beta-1)+(1:tn)).*VV(:,tn*(id_beta-1)+(1:tn))).*dN_int),2);
    end
    end
   sig_beta=diag(sqrt(pinv(A)*Sigma*pinv(A)));
   cov_beta= pinv(A)*Sigma*pinv(A) ;
   