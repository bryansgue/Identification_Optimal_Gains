function psid = adjust_angle(psid)
    while psid > 2*pi
        psid = psid - 4*pi;
    end
    while psid < -2*pi
        psid = psid + 4*pi;
    end
end
