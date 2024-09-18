% Expand a function and select all the code inside and use CTRL+T to
% uncomment. Only uncomment 1 at a time and run. Use CTRL+R to comment
% again after selecting all the code


%%%% Solutions for the 1st problem %%%%%%%

% function ExplicitEuler()
% 
%     clear all
%     close all
%     %%%Own Euler Solver
%     [t1,x1]=MyEuler(@dxdt,[0 4],[pi/18; 0],0.2);
%     
%     [t2,x2]=MyEuler(@dxdt,[0 4],[pi/18; 0],0.1);
%     
%     [t3,x3]=MyEuler(@dxdt,[0 4],[pi/18; 0],0.05);
%     
%     %%%Matlab Solver
%     [tmat,xmat]=ode23(@dxdt,[0 4],[pi/18; 0]);
%     
%     %%%plotting solution
%     figure(1)
%     plot(t1,x1(:,1),'k-',t2,x2(:,1),'b-',t3,x3(:,1),'m-',tmat,xmat(:,1),'ro-');
%     legend('h=0.2','h=0.1','h=0.05','ode23()');
%     xlabel('t');
%     ylabel('x');
%     figure(2)
%     plot(t1,x1(:,2),'k-',t2,x2(:,2),'b-',t3,x3(:,2),'m-',tmat,xmat(:,2),'ro-');
%     legend('h=0.2','h=0.1','h=0.05','ode23()');
%     xlabel('t');
%     ylabel('x');
%     
%     function [t,x]=MyEuler(ODEfunc,tspan,x0,h)
%         %%%This function solves system of ODEs.
%         t=tspan(1):h:tspan(2);
%         
%         x(1,:)=x0;
%         for n=1:length(t)-1
%             x(n+1,:)=x(n,:)+ODEfunc(t(n),x(n,:))'*h;
%         end
%     end 
%      
%     function xp=dxdt(tn,xn)
%         xp(1) = xn(2);
%         xp(2) = (-9.81/0.6)*xn(1);
%         xp = xp';
%     end
% end

% function ImplicitEuler()
%     clear all
%     
%     %%%%Solving using Implicit Euler
%     [t1,x1]=MyImpEuler(@dxdt,[0, 4],[pi/18, 0],0.2);
%     [t2,x2]=MyImpEuler(@dxdt,[0, 4],[pi/18, 0],0.1);
%     [t3,x3]=MyImpEuler(@dxdt,[0, 4],[pi/18, 0],0.05);
%     
%     figure(1)
%     plot(t1,x1(:,1),'k-',t2,x2(:,1),'b-',t3,x3(:,1),'m-');
%     legend('h=0.2','h=0.1','h=0.05','ode23()');
%     xlabel('t','Fontsize',14);
%     ylabel('x','Fontsize',14);
%     
%     figure(2)
%     plot(t1,x1(:,2),'k-',t2,x2(:,2),'b-',t3,x3(:,2),'m-');
%     legend('h=0.2','h=0.1','h=0.05','ode23()');
%     xlabel('t','Fontsize',14);
%     ylabel('x','Fontsize',14);
%     
%     function [t,x]=MyImpEuler(ODEfunc,tspan,x0,h)
%         t=tspan(1):h:tspan(2);
%         x(1,:)=x0;
%         
%         for n=1:length(t)-1
%             xn=x(n,:);
%             tnp1=t(n+1);
%             x(n+1,:)=fsolve(@(xnp1)EulerODEfunc(ODEfunc,tnp1,xn,xnp1,h),x(n,:));
%         end
%     end
%     
%     function residual=EulerODEfunc(ODE,tnp1,xn,xnp1,h)
%         residual=xn+h*ODE(tnp1,xnp1)-xnp1;
%     end
%     
%     function xp=dxdt(tn,xn)
%         xp(1) = xn(2);
%         xp(2) = (-9.81/0.6)*xn(1);
%     end
% 
% end

% function CrankNicolson()
%     clear all
%     
%     %%%%%Solving solution using different ODE
%     [t1,x1]=MyCrankNicolson(@dxdt,[0, 4],[pi/18, 0],0.2);
%     [t2,x2]=MyCrankNicolson(@dxdt,[0, 4],[pi/18, 0],0.1);
%     [t3,x3]=MyCrankNicolson(@dxdt,[0, 4],[pi/18, 0],0.05);
%      
%     %%%%Plotting solution
%     figure(1)
%     plot(t1,x1(:,1),'k-',t2,x2(:,1),'b-',t3,x3(:,1),'m-');
%     legend('h=0.2','h=0.1','h=0.05','ode23()');
%     xlabel('t','Fontsize',14);
%     ylabel('x','Fontsize',14);
%     
%     figure(2)
%     plot(t1,x1(:,2),'k-',t2,x2(:,2),'b-',t3,x3(:,2),'m-');
%     legend('h=0.2','h=0.1','h=0.05','ode23()');
%     xlabel('t','Fontsize',14);
%     ylabel('x','Fontsize',14);
%     
%     function [t,x]=MyCrankNicolson(ODEfunc,tspan,x0,h)
%         t=tspan(1):h:tspan(2);
%         x(1,:)=x0;
%         
%         for n=1:length(t)-1
%             xn=x(n,:);
%             tnp1=t(n+1);
%             tn=t(n);
%             x(n+1,:)=fsolve(@(xnp1)CrankNicolsonODEfunc(ODEfunc,tn,tnp1,xn,xnp1,h),x(n,:));
%         end
%     end
%     
%     function residual=CrankNicolsonODEfunc(ODE,tn,tnp1,xn,xnp1,h)
%         residual=xn+(h/2)*(ODE(tnp1,xnp1)+ODE(tn,xn))-xnp1;
%     end
%     
%     function xp=dxdt(tn,xn)
%         xp(1) = xn(2);
%         xp(2) = (-9.81/0.6)*xn(1);
%     end
% 
% end

% function RK2()
% 
%     close all
%     clear all
%     
%     %%%%%Own RK4 Solver
%     [t1,x1]=MyRK2sys(@dxdt,[0 4],[pi/18, 0],0.2);
%     [t2,x2]=MyRK2sys(@dxdt,[0 4],[pi/18, 0],0.1);
%     [t3,x3]=MyRK2sys(@dxdt,[0 4],[pi/18, 0],0.05);
%     %%%%%Matlab Solver
%     [tmat,xmat]=ode23(@dxdt,[0 4],[pi/18, 0]);
%     figure(1)
%     plot(t1,x1(:,1),'k-',t2,x2(:,1),'b-',t3,x3(:,1),'m-',tmat,xmat(:,1),'ro-');
%     legend('h=0.2','h=0.1','h=0.05','ode23()');
%     xlabel('t');
%     ylabel('x');
%     figure(2)
%     plot(t1,x1(:,2),'k-',t2,x2(:,2),'b-',t3,x3(:,2),'m-',tmat,xmat(:,2),'ro-');
%     legend('h=0.2','h=0.1','h=0.05','ode23()');
%     xlabel('t');
%     ylabel('x');
%     
%     function [t,x]=MyRK2sys(ODEfunc,tspan,x0,h)
%         %%%%%This function assumes you are only solving
%         %%%%one ODE.
%         t=tspan(1):h:tspan(2);
%         x(1,:)=x0;
%         for n=1:length(t)-1
%             k1=ODEfunc(t(n),x(n,:))';
%             k2=ODEfunc(t(n)+h,x(n,:)+h*k1)';
%             x(n+1,:)=x(n,:)+((h/2)*(k1+k2));
%         end
%     end
%      
%     function xp=dxdt(tn,xn)
%         xp(1) = xn(2);
%         xp(2) = (-9.81/0.6)*xn(1);
%         xp = xp';
%     end
% 
% end

% function RK4()
% 
%     close all
%     clear all
%     
%     %%%%Own RK4 Solver
%     [t1,x1]=MyRK4sys(@dxdt,[0 4],[pi/18, 0],0.2);
%     [t2,x2]=MyRK4sys(@dxdt,[0 4],[pi/18, 0],0.1);
%     [t3,x3]=MyRK4sys(@dxdt,[0 4],[pi/18, 0],0.05);
%     %%%%Matlab Solver
%     [tmat,xmat]=ode23(@dxdt,[0 4],[pi/18, 0]);
%     figure(1)
%     plot(t1,x1(:,1),'k-',t2,x2(:,1),'b-',t3,x3(:,1),'m-',tmat,xmat(:,1),'ro-');
%     legend('h=0.2','h=0.1','h=0.05','ode23()');
%     xlabel('t');
%     ylabel('x');
%     figure(2)
%     plot(t1,x1(:,2),'k-',t2,x2(:,2),'b-',t3,x3(:,2),'m-',tmat,xmat(:,2),'ro-');
%     legend('h=0.2','h=0.1','h=0.05','ode23()');
%     xlabel('t');
%     ylabel('x');
%     
%     
%     function [t,x]=MyRK4sys(ODEfunc,tspan,x0,h)
%         %%%%This function uses 4th order Runge-Kutta method to solve a
%         %%%%system of ODEs
%         t=tspan(1):h:tspan(2);
%         x(1,:)=x0;
%         for n=1:length(t)-1
%             k1=ODEfunc(t(n),x(n,:))';
%             k2=ODEfunc(t(n)+(h/2),x(n,:)+(h/2)*k1)';
%             k3=ODEfunc(t(n)+(h/2),x(n,:)+(h/2)*k2)';
%             k4=ODEfunc(t(n)+h,x(n,:)+h*k3)';
%             x(n+1,:)=x(n,:)+((1/6)*k1+(1/3)*(k2+k3)+(1/6)*k4)*h;
%         end
%     end
%      
%     function xp=dxdt(tn,xn)
%         xp(1) = xn(2);
%         xp(2) = (-9.81/0.6)*xn(1);
%         xp = xp';
%     end
% 
% end


%%%%%%% Solutions for the 2nd problem %%%%%%%%

% function ExplicitEuler2nd()
%     clear all
%     close all
%     %%%%%Own Euler Solver
%     [t1,x1]=MyEuler(@dxdt,[0 4],[pi/18; 0],0.2);
%     
%     [t2,x2]=MyEuler(@dxdt,[0 4],[pi/18; 0],0.1);
%     
%     [t3,x3]=MyEuler(@dxdt,[0 4],[pi/18; 0],0.05);
%     
%     %%%%Matlab Solver
%     [tmat,xmat]=ode23(@dxdt,[0 4],[pi/18; 0]);
%     
%     %%%%%plotting solution
%     figure(1)
%     plot(t1,x1(:,1),'k-',t2,x2(:,1),'b-',t3,x3(:,1),'m-',tmat,xmat(:,1),'ro-');
%     legend('h=0.2','h=0.1','h=0.05','ode23()');
%     xlabel('t');
%     ylabel('x');
%     figure(2)
%     plot(t1,x1(:,2),'k-',t2,x2(:,2),'b-',t3,x3(:,2),'m-',tmat,xmat(:,2),'ro-');
%     legend('h=0.2','h=0.1','h=0.05','ode23()');
%     xlabel('t');
%     ylabel('x');
%     
%     function [t,x]=MyEuler(ODEfunc,tspan,x0,h)
%         %This function solves system of ODEs.
%         t=tspan(1):h:tspan(2);
%         
%         x(1,:)=x0;
%         for n=1:length(t)-1
%             x(n+1,:)=x(n,:)+ODEfunc(t(n),x(n,:))'*h;
%         end
%     end 
%      
%     function xp=dxdt(tn,xn)
%         xp(1) = xn(2);
%         xp(2) = (-(1/2)*xn(2))-((49/3)*xn(1));
%         xp = xp';
%     end
% 
% end

% function ImplicitEuler2nd()
%     clear all
%     
%     %%%%%Solving using Implicit Euler
%     [t1,x1]=MyImpEuler(@dxdt,[0, 4],[pi/18, 0],0.2);
%     [t2,x2]=MyImpEuler(@dxdt,[0, 4],[pi/18, 0],0.1);
%     [t3,x3]=MyImpEuler(@dxdt,[0, 4],[pi/18, 0],0.05);
%     
%     figure(1)
%     plot(t1,x1(:,1),'k-',t2,x2(:,1),'b-',t3,x3(:,1),'m-');
%     legend('h=0.2','h=0.1','h=0.05','ode23()');
%     xlabel('t','Fontsize',14);
%     ylabel('x','Fontsize',14);
%     
%     figure(2)
%     plot(t1,x1(:,2),'k-',t2,x2(:,2),'b-',t3,x3(:,2),'m-');
%     legend('h=0.2','h=0.1','h=0.05','ode23()');
%     xlabel('t','Fontsize',14);
%     ylabel('x','Fontsize',14);
%     
%     function [t,x]=MyImpEuler(ODEfunc,tspan,x0,h)
%         t=tspan(1):h:tspan(2);
%         x(1,:)=x0;
%         
%         for n=1:length(t)-1
%             xn=x(n,:);
%             tnp1=t(n+1);
%             x(n+1,:)=fsolve(@(xnp1)EulerODEfunc(ODEfunc,tnp1,xn,xnp1,h),x(n,:));
%         end
%     end
%     
%     function residual=EulerODEfunc(ODE,tnp1,xn,xnp1,h)
%         residual=xn+h*ODE(tnp1,xnp1)-xnp1;
%     end
%     
%     function xp=dxdt(tn,xn)
%         xp(1) = xn(2);
%         xp(2) = (-(1/2)*xn(2))-((49/3)*xn(1));
%     end
% 
% end

% function CrankNicolson2nd()
%  
%     clear all
%     
%     %%%%%Solving solution using different ODE
%     [t1,x1]=MyCrankNicolson(@dxdt,[0, 4],[pi/18, 0],0.2);
%     [t2,x2]=MyCrankNicolson(@dxdt,[0, 4],[pi/18, 0],0.1);
%     [t3,x3]=MyCrankNicolson(@dxdt,[0, 4],[pi/18, 0],0.05);
%      
%     %%%%%Plotting solution
%     figure(1)
%     plot(t1,x1(:,1),'k-',t2,x2(:,1),'b-',t3,x3(:,1),'m-');
%     legend('h=0.2','h=0.1','h=0.05','ode23()');
%     xlabel('t','Fontsize',14);
%     ylabel('x','Fontsize',14);
%     
%     figure(2)
%     plot(t1,x1(:,2),'k-',t2,x2(:,2),'b-',t3,x3(:,2),'m-');
%     legend('h=0.2','h=0.1','h=0.05','ode23()');
%     xlabel('t','Fontsize',14);
%     ylabel('x','Fontsize',14);
%     
%     function [t,x]=MyCrankNicolson(ODEfunc,tspan,x0,h)
%         t=tspan(1):h:tspan(2);
%         x(1,:)=x0;
%         
%         for n=1:length(t)-1
%             xn=x(n,:);
%             tnp1=t(n+1);
%             tn=t(n);
%             x(n+1,:)=fsolve(@(xnp1)CrankNicolsonODEfunc(ODEfunc,tn,tnp1,xn,xnp1,h),x(n,:));
%         end
%     end
%     
%     function residual=CrankNicolsonODEfunc(ODE,tn,tnp1,xn,xnp1,h)
%         residual=xn+(h/2)*(ODE(tnp1,xnp1)+ODE(tn,xn))-xnp1;
%     end
%     
%     function xp=dxdt(tn,xn)
%         xp(1) = xn(2);
%         xp(2) = (-(1/2)*xn(2))-((49/3)*xn(1));
%     end
% 
% end

% function RK2_2nd()
%     
%     close all
%     clear all
%     
%     %%%%Own RK4 Solver
%     [t1,x1]=MyRK2sys(@dxdt,[0 4],[pi/18, 0],0.2);
%     [t2,x2]=MyRK2sys(@dxdt,[0 4],[pi/18, 0],0.1);
%     [t3,x3]=MyRK2sys(@dxdt,[0 4],[pi/18, 0],0.05);
%     %%%%Matlab Solver
%     [tmat,xmat]=ode23(@dxdt,[0 4],[pi/18, 0]);
%     figure(1)
%     plot(t1,x1(:,1),'k-',t2,x2(:,1),'b-',t3,x3(:,1),'m-',tmat,xmat(:,1),'ro-');
%     legend('h=0.2','h=0.1','h=0.05','ode23()');
%     xlabel('t');
%     ylabel('x');
%     figure(2)
%     plot(t1,x1(:,2),'k-',t2,x2(:,2),'b-',t3,x3(:,2),'m-',tmat,xmat(:,2),'ro-');
%     legend('h=0.2','h=0.1','h=0.05','ode23()');
%     xlabel('t');
%     ylabel('x');
%     
%     function [t,x]=MyRK2sys(ODEfunc,tspan,x0,h)
%         %%%%This function assumes you are only solving
%         %%%%one ODE.
%         t=tspan(1):h:tspan(2);
%         x(1,:)=x0;
%         for n=1:length(t)-1
%             k1=ODEfunc(t(n),x(n,:))';
%             k2=ODEfunc(t(n)+h,x(n,:)+h*k1)';
%             x(n+1,:)=x(n,:)+((h/2)*(k1+k2));
%         end
%     end
%      
%     function xp=dxdt(tn,xn)
%         xp(1) = xn(2);
%         xp(2) = (-(1/2)*xn(2))-((49/3)*xn(1));
%         xp = xp';
%     end
% 
% end

% function RK4_2nd()
%     
%     close all
%     clear all
%     
%     %%%%Own RK4 Solver
%     [t1,x1]=MyRK4sys(@dxdt,[0 4],[pi/18, 0],0.2);
%     [t2,x2]=MyRK4sys(@dxdt,[0 4],[pi/18, 0],0.1);
%     [t3,x3]=MyRK4sys(@dxdt,[0 4],[pi/18, 0],0.05);
%     %%%%%Matlab Solver
%     [tmat,xmat]=ode23(@dxdt,[0 4],[pi/18, 0]);
%     figure(1)
%     plot(t1,x1(:,1),'k-',t2,x2(:,1),'b-',t3,x3(:,1),'m-',tmat,xmat(:,1),'ro-');
%     legend('h=0.2','h=0.1','h=0.05','ode23()');
%     xlabel('t');
%     ylabel('x');
%     figure(2)
%     plot(t1,x1(:,2),'k-',t2,x2(:,2),'b-',t3,x3(:,2),'m-',tmat,xmat(:,2),'ro-');
%     legend('h=0.2','h=0.1','h=0.05','ode23()');
%     xlabel('t');
%     ylabel('x');
%     
%     
%     function [t,x]=MyRK4sys(ODEfunc,tspan,x0,h)
%         %%%%This function uses 4th order Runge-Kutta method to solve a
%         %%%%system of ODEs
%         t=tspan(1):h:tspan(2);
%         x(1,:)=x0;
%         for n=1:length(t)-1
%             k1=ODEfunc(t(n),x(n,:))';
%             k2=ODEfunc(t(n)+(h/2),x(n,:)+(h/2)*k1)';
%             k3=ODEfunc(t(n)+(h/2),x(n,:)+(h/2)*k2)';
%             k4=ODEfunc(t(n)+h,x(n,:)+h*k3)';
%             x(n+1,:)=x(n,:)+((1/6)*k1+(1/3)*(k2+k3)+(1/6)*k4)*h;
%         end
%     end
%      
%     function xp=dxdt(tn,xn)
%         xp(1) = xn(2);
%         xp(2) = (-(1/2)*xn(2))-((49/3)*xn(1));
%         xp = xp';
%     end
% 
% end