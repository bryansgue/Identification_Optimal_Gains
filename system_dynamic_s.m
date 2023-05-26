function [v] = system_dynamic_s(x, v, vref, h, L, ts)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

vp = dynamic_func(v,x , vref, L);
v = v +vp*ts;
end

