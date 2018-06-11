function [resi,predY]=residual(n,tn,alpha_n,beta_n,theta_n,gamma_n,alpha_hatT,beta_hat_wei,X,T,U,TU,VV ,Y,link_id,dN  )

    temp=0;
    for id_alpha=1:alpha_n
       temp= temp+alpha_hatT(:,tn*(id_alpha-1)+(1:tn)).*X(:,tn*(id_alpha-1)+(1:tn));
    end
    for id_beta=1:beta_n
	       zeta001(1:n,tn*(id_beta-1)+(1:tn))=repmat(beta_hat_wei(id_beta),n,tn);
	end
          
	[zeta002,dzeta002]=gamma_func(beta_hat_wei(beta_n+1:beta_n+theta_n),T,U,TU,theta_n,gamma_n);
	if beta_n>0  zeta00=[zeta001,zeta002]; else zeta00=[zeta002];end
    temp02=0;
    for id_beta=1:beta_n+theta_n;
	    temp02= temp02+zeta00(:,tn*(id_beta-1)+(1:tn)).*VV(:,tn*(id_beta-1)+(1:tn));
    end
    [predY,mu_dot_beta]=psy(temp+temp02,link_id);
    resi=(Y-predY).*dN;
