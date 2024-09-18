function Shooting_linear()

close all
clc

% Solving the ODE for
% 0 < x < 10 with x1(0)= 40 and x2(0) = z0 (Initial guess)
x0 = 0; xL = 10; h = 2; T0 = 40; TL = 200;

%Own RK4 Solver
z01 = 8;
[x1,T1]=MyRK4Sys(@dTdx,[x0 xL],[T0 z01],h)
disp(T1(end,1));

z02 = 10;
[x2,T2]=MyRK4Sys(@dTdx,[x0 xL],[T0 z02],h);
disp(T2(end,1));


z03 = z01 + (z02-z01)*(TL-T1(end,1))/(T2(end,1)-T1(end,1))
[x3,T3]=MyRK4Sys(@dTdx,[x0 xL],[T0 z03],h);
disp(T3(:,1));
%Matlab Solver
plot(x1,T1(:,1),'k-',x2,T2(:,1),'b-',x3,T3(:,1),'ro-')
legend('z0=8','z0=10','z0=opt');
xlabel('x');
ylabel('T');
end

function [x,T]=MyRK4Sys(ODEfunc,xspan,T0,h)
%This function uses 4th order Runge-Kutta method to solve a
%system of ODEs
x=xspan(1):h:xspan(2);
T(1,:)=T0;
for n=1:length(x)-1
    k1=ODEfunc(x(n),T(n,:))'
    k2=ODEfunc(x(n)+(h/2),T(n,:)+(h/2)*k1)'
    k3=ODEfunc(x(n)+(h/2),T(n,:)+(h/2)*k2)'
    k4=ODEfunc(x(n)+h,T(n,:)+h*k3)'
    T(n+1,:)=T(n,:)+((1/6)*k1+(1/3)*(k2+k3)+(1/6)*k4)*h
end
end
 
 
function Tp=dTdx(x,T)

Tp(1)=T(2);
Tp(2)=-0.01*(20-T(1));
Tp=Tp';
end