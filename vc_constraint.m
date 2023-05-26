function [c, ceq] = vc_constraint(u, vc_min, vc_max)
    c = [u(:) - vc_max; vc_min - u];
    ceq = [];
end