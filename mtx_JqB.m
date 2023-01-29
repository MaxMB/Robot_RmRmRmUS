function [JqB] = mtx_JqB (q4, q5, lB)
JqB = [(-1).*lB.*cos(q4).*cos(q5), lB.*sin(q4).*sin(q5);
    (-1).*lB.*cos(q5).*sin(q4), (-1).*lB.*cos(q4).*sin(q5);
    0, lB.*cos(q5)];