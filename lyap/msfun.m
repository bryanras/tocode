% this function defines the vector field of the DE to be solved
% as well as the linearized vector field and the tau-derivative
% differential equation
% The syntax is the following
% function [xp,Yp,vp]=msfun(tau,x,Y,v,f,Df)
%
% Input: 
% tau: value of period (current value)
% x: current value of solution
% Y: current value of transition matrix
% v: current value of tau-derivative
% f: string specifying the DE we want to solve
% Df: string specifying the Jacobian of the DE to solve
% Output:
% xp: derivative of x, vector field
% Yp: derivative for transition matrix
% vp: derivative for tau-derivative

function [xp,Yp,vp]=msfun(tau,x,Y,v,f,Df)
xp=feval(f,x); Jac=feval(Df,x); Yp=tau*Jac*Y; vp=tau*Jac*v+xp; xp=tau*xp; 
