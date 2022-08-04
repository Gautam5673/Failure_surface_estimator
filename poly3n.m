function [zz] = poly3n(param, x)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    x_max=param.f;
    x_min=param.g;
    zz=param.a+ param.e+ param.b*((x-x_min)./(x_max-x_min))+ param.c*((x-x_min)./(x_max-x_min)).^2+param.d*((x-x_min)./(x_max-x_min)).^3;

end