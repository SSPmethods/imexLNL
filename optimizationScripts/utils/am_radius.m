function r = am_radius(A,b)
%function r = am_radius(A,b)
%
%By David Ketcheson
%
%Evaluates the Radius of absolute monotonicity
%of a Runge-Kutta method, given the Butcher array.
%
%For an m-stage method, A should be an m x m matrix
%and b should be a column vector of length m.
%
%Accuracy can be changed by modifying the value of eps.
%Methods with very large radii of a.m. (>50) will require
%rmax to be increased.

rmax=50; eps=1.e-12;

m=length(b); e=ones(m,1);
K=[A;b'];
rlo=0; rhi=rmax;

while rhi-rlo>eps  %use bisection
  r=0.5*(rhi+rlo);
  X=eye(m)+r*A; beta=K/X; ech=r*K*(X\e);
  if (min(beta(:))<-3.e-16 || max(ech(:))>1.+3.e-16)
    rhi=r;
  else
    rlo=r;
  end
end

if rhi==rmax % r>=rmax
  error('Error: increase value of rmax in am_radius.m');
else
  r=rlo;
end
