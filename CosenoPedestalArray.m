function output = CosenoPedestalArray(SSL,p,N)

delta=0;
k=1;

Amplitudes = [];
delta1= (-13 + sqrt(13^2-(4*22*(13+SSL))))/44;
delta2= (-13 - sqrt(13^2-(4*22*(13+SSL))))/44;

if(delta1<0 && delta2<0)
    error('Ambas soluciones de la ecuación son negativas');
end
if(delta1<0)
    delta=delta2;
end
if(delta2<0)
    delta=delta1;
end

for n=1 : N
    Amplitud= delta + (1-delta)*(cos(((2*n-N-1)/(N-1))*pi/2))^p;
    Amplitudes(n)=Amplitud;
    k=k+1;
end
%plot(Amplitudes)
output = [Amplitudes(1);Amplitudes(2);Amplitudes(3);Amplitudes(4);Amplitudes(5);Amplitudes(6)];
end





