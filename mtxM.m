function [M] = mtxM (q1, q2, q3, q4, q5, lG2, l3, lG3, l4, J1, J2, ...
    J3, m2, m3)
M = [J1+lG2.^2.*m2.*sin(q2).^2+m3.*(l3.*sin(q2)+lG3.*sin(q2+q3)).^2+ ...
  J3.*l4.^(-2).*cos(q5).^2.*(l3.*sin(q2)+l4.*sin(q2+q3)).^2.*sin(q1+ ...
  (-1).*q4).^2.*(cos(q1).*cos(q2+q3).*cos(q4).*cos(q5)+cos(q2+q3).* ...
  cos(q5).*sin(q1).*sin(q4)+sin(q2+q3).*sin(q5)).^(-2)+l4.^(-2).* ...
  lG3.^2.*m3.*cos(q5).^2.*(l3.*sin(q2)+l4.*sin(q2+q3)).^2.*sin(q1+( ...
  -1).*q4).^2.*(cos(q1).*cos(q2+q3).*cos(q4).*cos(q5)+cos(q2+q3).* ...
  cos(q5).*sin(q1).*sin(q4)+sin(q2+q3).*sin(q5)).^(-2), ...
  l3.*l4.^(-2) ...
  .*cos(q5).*(l3.*sin(q2)+l4.*sin(q2+q3)).*sin(q1+(-1).*q4).*(cos( ...
  q1).*cos(q2+q3).*cos(q4).*cos(q5)+cos(q2+q3).*cos(q5).*sin(q1).* ...
  sin(q4)+sin(q2+q3).*sin(q5)).^(-2).*((-1).*J3.*(cos(q1).*cos(q2).* ...
  cos(q4).*cos(q5)+cos(q2).*cos(q5).*sin(q1).*sin(q4)+sin(q2).*sin( ...
  q5))+lG3.*m3.*((-1).*lG3.*(cos(q1).*cos(q2).*cos(q4).*cos(q5)+cos( ...
  q2).*cos(q5).*sin(q1).*sin(q4)+sin(q2).*sin(q5))+l4.*cos(q3).*( ...
  cos(q1).*cos(q2+q3).*cos(q4).*cos(q5)+cos(q2+q3).*cos(q5).*sin(q1) ...
  .*sin(q4)+sin(q2+q3).*sin(q5))));
  l3.*l4.^(-2).*cos(q5).*(l3.*sin(q2)+ ...
  l4.*sin(q2+q3)).*sin(q1+(-1).*q4).*(cos(q1).*cos(q2+q3).*cos( ...
  q4).*cos(q5)+cos(q2+q3).*cos(q5).*sin(q1).*sin(q4)+sin(q2+q3).* ...
  sin(q5)).^(-2).*((-1).*J3.*(cos(q1).*cos(q2).*cos(q4).*cos(q5)+ ...
  cos(q2).*cos(q5).*sin(q1).*sin(q4)+sin(q2).*sin(q5))+lG3.*m3.*(( ...
  -1).*lG3.*(cos(q1).*cos(q2).*cos(q4).*cos(q5)+cos(q2).*cos(q5).* ...
  sin(q1).*sin(q4)+sin(q2).*sin(q5))+l4.*cos(q3).*(cos(q1).*cos(q2+ ...
  q3).*cos(q4).*cos(q5)+cos(q2+q3).*cos(q5).*sin(q1).*sin(q4)+sin( ...
  q2+q3).*sin(q5)))), ...
  J2+lG2.^2.*m2+l3.^2.*m3.*sin(q3).^2+J3.*l3.^2.* ...
  l4.^(-2).*(cos(q1).*cos(q2).*cos(q4).*cos(q5)+cos(q2).*cos(q5).* ...
  sin(q1).*sin(q4)+sin(q2).*sin(q5)).^2.*(cos(q1).*cos(q2+q3).*cos( ...
  q4).*cos(q5)+cos(q2+q3).*cos(q5).*sin(q1).*sin(q4)+sin(q2+q3).* ...
  sin(q5)).^(-2)+l3.^2.*m3.*(cos(q3)+(-1).*l4.^(-1).*lG3.*(cos(q1).* ...
  cos(q2).*cos(q4).*cos(q5)+cos(q2).*cos(q5).*sin(q1).*sin(q4)+sin( ...
  q2).*sin(q5)).*(cos(q1).*cos(q2+q3).*cos(q4).*cos(q5)+cos(q2+q3).* ...
  cos(q5).*sin(q1).*sin(q4)+sin(q2+q3).*sin(q5)).^(-1)).^2];