clear;
clc;

% ANALISIS DE LA BOCINA 
%Parametros de entrada
A = 0.21; %Anchura de la bocina
B = 0.76; %Altura de la bocina
s = 1/4; %Error de fase plano E

f = 1.060; %Frecuencia de trabajo

%Calculos previos
l_onda = 0.3/f; %Calculo de la longitud de onda
X =[]; %Array donde almacenamos los distintos valores de anchura de la bocina
Y =[]; %Array donde almacenamos los distintos valores de altura de la bocina
E_O=1; 
i_ancho_x = 1; %Indice que nos permite recorrer el ancho de la bocina en el bucle
i_alto_y = 1; %Indice que nos permite recorrer el alto de la bocina en el bucle
saltos=100; %muestreo la superficie

%Calculo de R1 Y R2 utiles para el campo de la bocina
%R1 = (A^2)/(8*l_onda*t)
R2 = (B^2)/(8*l_onda*s)

%Array que almacena los distintos valores del ancho de la bocina
for x = linspace(-A/2,A/2,saltos)  
    X(i_ancho_x) = x;
    i_ancho_x=i_ancho_x+1;
end

%Array que almacena los distintos valores del alto de la bocina
for y = linspace(-B/2,B/2,saltos) 
    Y(i_alto_y)=y;
    i_alto_y = i_alto_y+1;

end

%Campo de apertura de la bocina
E_y_mod = E_O*cos((pi*X)/A); %Modulo del campo
beta = (2*pi/l_onda)*sqrt(1-(l_onda/(2*A))^2);    %Variable util para calcular la fase 
E_y_fase = -beta*1i*((Y.^2)/(2*R2));    %Fase del campo
E_y_tot = E_y_mod.*exp(E_y_fase); %

%Campo total de apertura de la bocina

% Cortes de la bocina en phi=0 y phi=pi/2
theta = [];
theta_rad=[];
Py_bocina = []; %factor de radiacion de la bocina
Py_sector_bocina =[]; %
Py_sector_con_factor_Array=[];
aux = 0; %valor maximo en la maxima direccion, para normalizar

%Ancho de haz en el eje horizontal, dimension A 
PHI = 0; %El eje horizontal cuenta con phi=0
i_THETA = 1; %Indice que nos permitira recorre los distintos valores de theta 

%Recorremos los distintos valores de theta para phi=0
 for THETA = linspace(-pi/2,pi/2,500+1)
       campo_ap=@(x,y) E_O.*cos((pi.*x)./A).*exp(-beta.*1i.*(((y.^2)./(2.*R2)))).*exp(1i*(2.*pi./l_onda).*(sin(THETA).*cos(PHI).*x+sin(THETA).*sin(PHI).*y));
       E_rad_boc = integral2(campo_ap,-A/2,A/2,-B/2,B/2);
       
       if (THETA == 0) %Direccion de maxima radiaci?n de la bocina
       aux = abs(E_rad_boc); 
       end
          
       Py_bocina(i_THETA)= E_rad_boc; %Factor de radiacion de la bocina
       THETA_grados = THETA*180/pi;
       theta(i_THETA) = THETA_grados;
       i_THETA =i_THETA+1;
 end
 
 FACTOR_ARRAY = AFCU_DIFERENCIA(1.06e9);
 CP= CosenoPedestalArray(-15,2,6);
 
 
 E_campo_diferencia = Py_bocina.* (CP(1,:).* FACTOR_ARRAY(1,:) + ...
                                   CP(2,:).* FACTOR_ARRAY(2,:) + ...
                                   CP(3,:).* FACTOR_ARRAY(3,:) + ...
            exp(j*pi).*(           CP(4,:).* FACTOR_ARRAY(4,:) + ...
                                   CP(5,:).* FACTOR_ARRAY(5,:) + ...
                                   CP(6,:).* FACTOR_ARRAY(6,:)));
 
                               
 E_campo_sector_horizontal =  Py_bocina.*AFCU(1.06e9);               
 E_campo_sector_horizontal_abs = abs(E_campo_sector_horizontal);
 Max_campo_horizontal = max(E_campo_sector_horizontal_abs);
 
 E_campo_diferencia_abs = abs(sum(E_campo_diferencia,1));
 Max_campo_diferencia = max(abs(E_campo_diferencia_abs));

 E_campo_suma = Py_bocina.*( FACTOR_ARRAY(1,:)+ FACTOR_ARRAY(2,:) + FACTOR_ARRAY(3,:) + (FACTOR_ARRAY(4,:) + FACTOR_ARRAY(5,:) + FACTOR_ARRAY(6,:)));
 E_campo_suma_abs = abs(sum(E_campo_suma,1));
 Max_campo_suma = max(abs(E_campo_suma_abs));

Py_mod= []; %Factor de radiacion normalizado para phi=0
Max_diagrama_diferencia =[];
Max_diagrama_suma=[];
Max_diagrama_sector =[];

i_ARRAY = 0; %Indice para recorrer el array

%Normalizamos el modulo del factor de radiaci?n
for ARRAY = linspace(1,500+1,500+1)
    
    norm = abs(Py_bocina(ARRAY))/aux;
    norm3 = abs(E_campo_diferencia_abs(ARRAY))/Max_campo_diferencia;
    norm4 = abs(E_campo_suma_abs(ARRAY))/Max_campo_suma;
    norm5 = abs(E_campo_sector_horizontal_abs(ARRAY))/Max_campo_horizontal;
    
    Py_mod(ARRAY) = 20*log(norm); %Factor de radiacion normalizado de la bocina
    Max_diagrama_diferencia(ARRAY) = 20*log(norm3);
    Max_diagrama_suma(ARRAY)= 20*log(norm4);
    Max_diagrama_sector(ARRAY)= 20*log(norm5);
    
    i_ARRAY = i_ARRAY+1;
end

%Graficas para el corte phi=0
figure(1)
set(gcf,'Name','Corte del plano horizontal de la bocina')
subplot(121)
plot(theta,Py_mod);
grid on
title('Factor de radiacion normalizado para \phi =0')
xlabel('\theta (grados)')
ylabel(' dB')
subplot(122)
theta_rad = linspace(-pi/2,pi/2,500+1);
polar(theta_rad,abs(Py_bocina/aux))
title('diagrama de radiacion normalizado para \phi =0')
grid on


%Grafica del diagrama del diferencia.
figure(2)
set(gcf,'Name','Grafica del diagrama del campo de diferencia')
subplot(121)
plot(theta,Max_diagrama_diferencia);
axis([-100 100 -40 0]);
grid on
title('Grafica del diagrama del campo de diferencia normalizado del un sector para \phi =0')
xlabel('\theta (grados)')
ylabel('diagrama diferencia mod(dB)')
subplot(122)
theta_rad = linspace(-pi/2,pi/2,500+1);
polar(theta_rad,E_campo_diferencia_abs)
title('Diagrama del diferencia \phi =0')
grid on


%Grafica del diagrama del suma.
figure(3)
set(gcf,'Name','Grafica del diagrama del campo de suma')
subplot(121)
plot(theta,Max_diagrama_suma);
%axis([-100 100 -40 0]);
grid on
title('Factor de radiacion normalizado del un sector para \phi =0')
xlabel('\theta (grados)')
ylabel('diagrama suma mod(dB)')
subplot(122)
theta_rad = linspace(-pi/2,pi/2,500+1);
polar(theta_rad,E_campo_suma_abs)
title('Diagrama del suma \phi =0')
grid on

%Grafica del diagrama del suma y diferencia.
figure(4)
set(gcf,'Name','Grafica del diagrama  suma y diferencia')
subplot(121)
plot(theta,Max_diagrama_suma);
hold on;
plot(theta,Max_diagrama_diferencia);
axis([-100 100 -40 0]);
grid on
title('diagrama del suma y resta del un sector para \phi =0')
xlabel('\theta (grados)')
ylabel('diagrama suma y resta (dB)')
subplot(122)
theta_rad = linspace(-pi/2,pi/2,500+1);
polar(theta_rad,E_campo_suma_abs)
hold on;
polar(theta_rad,E_campo_diferencia_abs)
title('Diagrama del suma \phi =0')
grid on

%Grafica del diagrama horizontal sector.
figure(5)
set(gcf,'Name','Grafica en el plano horizontal de un sector')
subplot(121)
plot(theta,Max_diagrama_sector);
hold on;
grid on
title('diagrama en plano horizontal de un sector para \phi =0')
xlabel('\theta (grados)')
ylabel(' dB')
subplot(122)
theta_rad = linspace(-pi/2,pi/2,500+1);
polar(theta_rad,E_campo_sector_horizontal_abs)
hold on;
title('Diagrama en el plano vertical \phi =0')
grid on

