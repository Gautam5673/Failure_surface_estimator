function [ param ] = spline_param_new(spm)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    param.a=spm.Y_1;
    param.b=spm.D_1;
    param.c=3*(spm.Y_2-spm.Y_1)-2*spm.D_1-spm.D_2;
    param.d=-2*(spm.Y_2-spm.Y_1)+spm.D_1+spm.D_2;
    param.e=spm.shift;
    param.f=spm.X_2;
    param.g=spm.X_1;
end
