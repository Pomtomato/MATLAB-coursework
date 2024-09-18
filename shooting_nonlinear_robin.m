function shooting_nonlinear_robin()

close all; clc;
init_guess = .005;
x1 = fzero(@solver, init_guess);
analyticalSOL();

    function F = solver(y01)
        x0 = 0; xL = 1; h = 0.5;

%         options=odeset('RelTol',1e-08,'AbsTol',[1e-08, 1e-08]);
%         [x, ode] = ode45(@myfunction,[x0, xL],[y01, -.5+y01], options);
%         
%         F = 1 - (ode(end,1) + ode(end,2)); %residue
%         figure(1)
%         plot(x,ode(:,1))
%         figure(2)
%         plot(x,ode(:,2))
        
        [x, RK] = MyRK4Sys(@myfunction, [x0 xL], [y01, -.5+y01], h);

        F = 1 - (RK(end,1) + RK(end,2)) ;%residue
        figure(1);
        plot(x,RK(:,1));
        figure(2);
        plot(x,RK(:,2));

    end

    function [x,Y]=MyRK4Sys(ODEfunc,xspan,Y0,h)

        x=xspan(1):h:xspan(2);
        Y(1,:)=Y0;
        for n=1:length(x)-1
            k1=ODEfunc(x(n),Y(n,:));
            k2=ODEfunc(x(n)+(h/2),Y(n,:)+(h/2)*k1);
            k3=ODEfunc(x(n)+(h/2),Y(n,:)+(h/2)*k2);
            k4=ODEfunc(x(n)+h,Y(n,:)+h*k3);
            Y(n+1,:)=Y(n,:)+((1/6)*k1+(1/3)*(k2+k3)+(1/6)*k4)*h;
        end
    end

    function dy = myfunction(x,Z)
        dy = zeros(2,1);
        dy(1) = Z(2);
        dy(2) = 0.5*(1+x+Z(1)).^3;
        dy = dy';
    end

    function analyticalSOL()

        x = linspace(0,1,3);
        y = zeros(size(x));
        for n = 1:numel(x)
          y(n) = (2/(2-x(n)))-x(n)-1;
        end
        figure(3);
        plot(x,y);
    end

end