%k-method - V-g method Flutter Calculation
%stead-flow aerodynamic solution

%V-g method flutter speed prediction
% Incorrectness of the k Method for Flutter Calculations 
% ZAERO AE Theory

close all;clear all;clc

%%
%Actual parameters of the system

b = 2; 
I_p = 1E+6;
m = 10;
k_t = 2e7;
k_h = 3e6;
rho = 2.23;

%%
%Non-dimensional parameters
e = 0.25;
a = 0;
sigma = sqrt(2);
r = sqrt(1/3);
x_theta = 0.25;
nu = 200;
w_h = sqrt(k_h/m);
%%
for U=1:100000
    

    V = U./(b*w_h);

    mat = inv([1 x_theta;x_theta r^2]) * [-1 -2*V^2/nu; 0 -sigma^2*r^2+2*V^2/nu*(0.5+a)];
    sol = eig(mat);
    fre(U) = sol(1,1);
    damp(U) = sol(2,1);
    Vel(U) = V;
    
end

figure(1)
plot(Vel,fre)
title("Velocity-Frequency")

figure(2)
plot(Vel,damp)
title("Velocity-Damping")

figure(3)
plot(Vel,fre,Vel,damp)
title("k method")
legend("Frequency","Damping")

