function [c, ceq] = h_p_constraint(h_p, h_p_min, h_p_max)
    % Restricción de variación temporal de hd_p
    delta_h_p = diff(h_p);
    c = [delta_h_p - h_p_max; -delta_h_p + h_p_min];
    ceq = [];
end