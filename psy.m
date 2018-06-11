function [val,div]=psy(x,link_id)
%Important!!  
%if this file is changed, the corresponding function in Pgamma_b_X.c must
%be changed in the same way,and it must be compiled again.



if link_id==0  %%%psy(x)=exp(x);
    val=exp(x);
    div=exp(x);
end

if link_id==1  %%%psy(x)=x;
    val=x;
    div=x-x+1;
end

if link_id==2  %%%psy(x)=logistic;
    val=exp(x)./(1+exp(x));
    div=exp(x)./((1+exp(x)).*(1+exp(x)));
end