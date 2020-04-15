clc; clear all; close all;

x0 = zeros(1,36);
fun = @areklmn_reactor;
x = fsolve(fun, [x0, 0, 0])% or x = solve(@func_example, x0), x0 are the initial guess vector, x is the outpur solution vector
% Extracting elements from x to get
