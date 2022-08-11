% Function to find F_lambda (used in pseudo arclength)
% Partial derivative of F wrt lambda

% Inputs:
% x - vector of x ordinates
% gamma - parameter (not used)

% Ouputs:
% F_lambda - vector of partial derivative of F wrt lambda

function F_lambda=findF_lambda(x,lambda)

F_lambda=sech(x).^2;

end