function [AF] = Arraycircularplano(f)

c=3e8;
l=c/f;
k=(2*pi)/l;
N=6;
n=0:N-1;
phi_n=2*pi*n/N;
phi=(0*pi)/180:((2*pi))/500:(2*pi);
M=length(phi);
d_circular=0.7*l;
a=(N*d_circular)/(2*pi);

for m=1:M

    AF(m)=sum(exp(i*k*a*(cos(phi(m)-phi_n))));
    
  
    
end
   

polar(phi, abs(AF))




end

