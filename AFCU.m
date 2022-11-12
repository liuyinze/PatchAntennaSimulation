function [AF] = AFCU(f)

c=3e8;
l=c/f;
k=(2*pi)/l;
N=6;
n=0:N-1;
phi_n=0:12*pi/180:72*pi/180;
M = length(phi_n);
d_circular=0.7*l;
r= 1.26;
phi = [0:pi/3:2*pi];
theta = [-pi/2:pi/500:pi/2];
THETA=  [-pi/2:pi/500:pi/2];
theta_grados=[];


    AFUCA = zeros(N,length(theta))
    for m =1:N
        for i_theta = 1:length(theta)
            AFUCA(m,i_theta) = sum(exp(j*k*r*( sin(theta(i_theta))*cos(0-phi_n(m)) )));
    
              THETA_grados = theta(i_theta)*180/pi;
       theta_grados(i_theta) = THETA_grados; 
        end
   
    end
    
AFUCA_ABS = abs(sum(AFUCA,1));
AF=sum(AFUCA,1);
MAX =  max(AF);
% plot(theta_grados,20*log(AFUCA_ABS/MAX));
% axis([-60 60 -80 0]);
% title('Fctor del array en el plano vertical');
% xlabel('\theta (grados)');
% ylabel('dB');
% grid on;


   




end




