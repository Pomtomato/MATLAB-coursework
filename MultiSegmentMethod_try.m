A = [0 1;-1 0];
B = [0;-1];
C =[1 0;0 0];
D = [0 0; 1 0];
E = [0 ;0];
multiSegmentEx1(A,B,C,D,E)

function multiSegmentEx1(A,B,C,D,E)
clc;
h = 0.5;

% A = [0 1;-1 0];
% B = [0;-1];
% C =[1 0;0 0];
% D = [0 0; 1 0];
% E = [0 ;0];

Y0 = [1 0;0 1];
Z0 = [0; 0];
x0 = 0; xL=1;

[x1,Y1]=RK4solver(@dYdx,[x0 xL],Y0(:,1),h,A,B);
[x1,Y2]=RK4solver(@dYdx,[x0 xL],Y0(:,2),h,A,B);
[x1,Z]=RK4solver(@dZdx,[x0 xL],Z0,h,A,B);

Ymid = [Y1(:,2) Y2(:,2)];
Yend = [Y1(:,3) Y2(:,3)];
Zmid = Z(:,2);
Zend = Z(:,3);

y_at_0 = inv(C + D*Yend)*(E - D*Zend);
y_at_half = Ymid * y_at_0 + Zmid;


function [x, Y] = RK4solver(myFunc, xspan, Y0, h,A,B)

x = xspan(1):h:xspan(2);
x = x';
%Ycol = length(Y0(1,:));
Y = Y0(:,1);
    for n = 1:length(x)-1
        k1 = myFunc(x(n), Y(:,n) ,A,B);
        k2 = myFunc(x(n)+(h/2), Y(:,n)+(h/2)*k1 ,A,B);
        k3 = myFunc(x(n)+(h/2), Y(:,n)+(h/2)*k2 ,A,B);
        k4 = myFunc(x(n)+h, Y(:,n)+h*k3 ,A,B);
        Y(:,n+1) = Y(:,n)+ ((1/6)*k1+(1/3)*(k2+k3)+(1/6)*k4)*h;
    end
end
% end

function dY = dYdx(x,Y,A,B)

n = 1;

dY(1) = A(1,:)*Y(:,n);
dY(2)= A(2,:)*Y(:,n);
dY = dY';
n= n+1;
end

function dZ = dZdx(x,Y,A,B)

n = 1;

dZ(1) = A(1,:)*Y(:,n) + B(1);
dZ(2)= A(2,:)*Y(:,n) + B(2);
dZ = dZ';
n= n+1;
end

end
