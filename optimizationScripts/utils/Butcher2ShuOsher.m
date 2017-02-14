function [Re,P,Q] = Butcher2ShuOsher(A,Ahat,b,bhat,r,r2);
%function [Re,P,Q] = Butcher2ShuOsher(A,Ahat,b,bhat,r,K);
% Re = v
% P = alpha
% Q = alpha_hat
%Converting Butcher form to Modified Shu Osher
%Y= Un+dt*A*F(Un)+dt^2*Ahat*Fdot(Un)
%Y= R*e*Un+ P*(Un+(dt/r)*F(Un))+Q*(Un+(dt^2/r2)*Fdot(Un))     
    
    %r2=r*K; 
    s=length(A); z=zeros(s+1,1); I=eye(s+1);e=ones(s+1,1);
	S=[[A;b],z];Shat=[[Ahat;bhat],z];          
    %Converting butcher to Modified Shu Osher
    %R=inv((I+r*S+r2*Shat));  
    
    Re=(I+r*S+r2*Shat)\e;
    P=(I+r*S+r2*Shat)\(r*S);
    Q=(I+r*S+r2*Shat)\(r2*Shat);
    
end
