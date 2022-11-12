function output = AFCU_DIFERENCIA(f)
c=3e8;
l=c/f;
k=(2*pi)/l;
N=6;
n=0:N-1;
phi_n=0:12*pi/180:72*pi/180;
M = length(phi_n);
d_circular=0.7*l;
r= 1.26;
theta = [-pi/2:pi/500:pi/2];
THETA=  [-pi/2:pi/500:pi/2];
theta_grado =[];

    AFUCA = zeros(N,length(theta));
    for m =1:N;
        for i_theta = 1:length(theta);
            AFUCA(m,i_theta) = sum(exp(j*k*r*( sin(theta(i_theta))*cos(0-phi_n(m)) )));
            theta_grado(i_theta)=theta(i_theta)*(180/pi);
        end
        
    end
    E1 = AFUCA(1,:);
    E2 = AFUCA(2,:);
    E3 = AFUCA(3,:);
    E4 = AFUCA(4,:);
    E5 = AFUCA(5,:);
    E6 = AFUCA(6,:);

    output = [E1;E2;E3;E4;E5;E6];
end

