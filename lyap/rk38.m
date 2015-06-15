% this function takes one integration step with RK-38 rule
% of the differential equations to be solved to set up
% the multiple shooting matrix
% The syntax is the following
% function [x1,Y1,v1]=rk38(h,tau,x0,Y0,v0,f,Df)
%
% Input: 
% h: stepsize
% tau: value of period (current value)
% x0: value of solution at beginning of step
% Y0: value of transition matrix at beginning of step
% v0: value of tau-derivative at beginning of step
% f: string specifying the DE we want to solve
% Df: string specifying the Jacobian of the DE to solve
% Output:
% x1: value of solution at end of step
% Y1: value of transition matrix at end of step
% v1: value of tau-derivative at end of step

function [x1,Y1,v1]=rk38(h,tau,x0,Y0,v0,f,Df)
[xk1,Yk1,vk1]=msfun(tau,x0,Y0,v0,f,Df); % 
xs=x0+(1/3)*h*xk1; Ys=Y0+(1/3)*h*Yk1;vs=v0+(1/3)*h*vk1; % stage values
[xk2,Yk2,vk2]=msfun(tau,xs,Ys,vs,f,Df); % 
xs=x0+h*(-(1/3)*xk1+xk2); Ys=Y0+h*(-(1/3)*Yk1+Yk2);vs=v0+h*(-(1/3)*vk1+vk2); % 
[xk3,Yk3,vk3]=msfun(tau,xs,Ys,vs,f,Df); % 
xs=x0+h*(xk1-xk2+xk3); Ys=Y0+h*(Yk1-Yk2+Yk3);vs=v0+h*(vk1-vk2+vk3); % 
[xk4,Yk4,vk4]=msfun(tau,xs,Ys,vs,f,Df); % 
x1=x0+h*(0.125*xk1+0.375*xk2+0.375*xk3+0.125*xk4);
Y1=Y0+h*(0.125*Yk1+0.375*Yk2+0.375*Yk3+0.125*Yk4);
v1=v0+h*(0.125*vk1+0.375*vk2+0.375*vk3+0.125*vk4);
