function [G] = vecG_RRR (q2, q3, lG2, l3, lG3, m2, m3, g)
G = [0;
    (-1).*g.*(lG2.*m2.*sin(q2)+m3.*(l3.*sin(q2)+lG3.*sin(q2+q3)));
    (-1).*g.*lG3.*m3.*sin(q2+q3)];
