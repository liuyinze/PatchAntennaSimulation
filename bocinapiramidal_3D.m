function [X,Y,Z,f] = bocinapiramidal_3D(a,b)
%BOCINA PIRAMIDAL 3D
% Definimos los �ngulos que definen cada direcci�n %
dt=pi/150; % Paso en theta %
dp=2*pi/150; % Paso en phi %
t=0:dt:pi/2; % �ngulo theta %
p=0:dp:2*pi; % �ngulo phi %
let=length(t); % N�mero de elementosde t %
lp=length(p); % N�mero de elementosde p %
k=2*pi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definimos los par�metros de la antena concreta %
%a=2; % Longitud el�ctricade a %
%b=2; % Longitud el�ctricade b %
% Definimos los campos de la antena concreta %
Et=zeros(let,lp); % Inicializamos a cerocomponente Et %
Ep=zeros(let,lp); % Inicializamos a cerocomponente Ep %
for I=1:let
 theta=t(I);
 for J=1:lp
 phi=p(J);
 u=(k*a/2).*sin(theta).*cos(phi);
 v=(k*b/2)*sin(theta)*sin(phi);

 d1=(pi/2)^2-u.^2;

 if (abs(d1)<=10^(-5))
 D1=1/pi;
 else
 D1=cos(u)/d1;
 end

 D2=sinc(v/pi);

 % Et(I,J)=(sin(phi)).*(sin(v)/v).*D1;
 % Ep(I,J)=(cos(phi)).*(cos(theta)).*(sin(v)/v).*D1;
 Et(I,J)=(sin(phi)).*D2.*D1;
 Ep(I,J)=(cos(phi)).*(cos(theta)).*D2.*D1;


 end
end
E=sqrt(abs(Et.^2)+abs(Ep.^2)); % M�dulo del campo electrico total %Desarrollo de una GUI para la representaci�n 2D y 3D del diagrama de radiaci�n de antenas P�gina 98
En=E/max(max(E)); % M�dulo del campoelectrico total NORMALIZADO %
Pn=En.^2; % Potencia normalizadaen escala lineal %
PndB=10*log10(Pn); % Potencia normalizada en dBs %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definimos el diagrama 3D %
X=zeros(let,lp); % Inicializamos X %
Y=zeros(let,lp); % Inicializamos Y %
Z=zeros(let,lp); % Inicializamos Z %
F=PndB; % Factor F %
Fmin=-40; % Valor m�nimo en dB %
f=0.5*(1-F/Fmin+abs(1-F/Fmin)); % Par�metro para eltrazado %
for I=1:let
 theta=t(I);
 for J=1:lp
 phi=p(J);
 X(I,J)=f(I,J)*sin(theta)*cos(phi);
 Y(I,J)=f(I,J)*sin(theta)*sin(phi);
 Z(I,J)=f(I,J)*cos(theta);
 end
 %surf(X,Y,Z,f)
 surf(X,Y,Z,f)
 %plot(E)
 %polar(theta,E)
end
