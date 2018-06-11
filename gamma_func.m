function [gamma,dgamma]=gamma_func(theta,T,U,TU,theta_n,gamma_n)
n=size(T,1);
tn=size(T,2);
%1st RANDONMIZATION
    %linear
         %  gamma=[repmat(theta(1),n,tn)+theta(4)*TU, repmat(theta(2),n,tn)+theta(5)*TU,repmat(theta(3),n,tn)+theta(6)*TU,repmat(0,n,tn*(theta_n-gamma_n))];
         %  dgamma=[repmat(1,n,tn),repmat(1,n,tn),repmat(1,n,tn), TU , TU , TU];
	   
	%quadratic	
         %  gamma=[repmat(theta(1),n,tn)+theta(4)*TU+theta(7)*TU.^2, repmat(theta(2),n,tn)+theta(5)*TU+theta(8)*TU.^2, repmat(theta(3),n,tn)+theta(6)*TU+theta(9)*TU.^2,repmat(0,n,tn*(theta_n-gamma_n))];
         %  dgamma=[repmat(1,n,tn),repmat(1,n,tn),repmat(1,n,tn), TU , TU , TU, TU.^2 , TU.^2 , TU.^2];
	
    %expdouble	
	     %	gamma=[repmat(theta(1),n,tn)+theta(4)*TU, repmat(theta(2),n,tn)+theta(5)*TU,repmat(theta(3),n,tn)+theta(6)*TU,repmat(0,n,tn*(theta_n-gamma_n))];
         % dgamma=[repmat(1,n,tn),repmat(1,n,tn),repmat(1,n,tn), TU , TU , TU];
	
   
   
%2nd RANDONMIZATION
         %   gamma=[repmat(theta(1),n,tn)+theta(3)*TU, repmat(theta(2),n,tn)+theta(4)*TU,repmat(0,n,tn*(theta_n-gamma_n))];
         %   dgamma=[repmat(1,n,tn),repmat(1,n,tn),TU , TU ];
	 
	%quadratic	
  
         gamma=[repmat(theta(1),n,tn)+theta(3)*TU+theta(5)*TU.*TU,repmat(theta(2),n,tn)+theta(4)*TU+theta(6)*TU.*TU ,repmat(0,n,tn*(theta_n-gamma_n))];
         dgamma=[repmat(1,n,tn),repmat(1,n,tn),TU,TU,TU.*TU,TU.*TU];
  

	


  
