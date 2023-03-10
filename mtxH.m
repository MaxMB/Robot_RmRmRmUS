function [H] = mtxH (q1, q2, q3, q4, q5, dq1, dq2, dq3, dq4, dq5, lG2, ...
    l3, lG3, l4, J3, m2, m3)
H = [dq2.*lG2.^2.*m2.*cos(q2).*sin(q2)+m3.*(l3.*sin(q2)+lG3.*sin(q2+ ...
  q3)).*(dq2.*l3.*cos(q2)+(dq2+dq3).*lG3.*cos(q2+q3)+dq1.*l4.^(-1).* ...
  lG3.*cos(q2+q3).*cos(q5).*(l3.*sin(q2)+l4.*sin(q2+q3)).*sin(q1+( ...
  -1).*q4).*(cos(q1).*cos(q2+q3).*cos(q4).*cos(q5)+cos(q2+q3).*cos( ...
  q5).*sin(q1).*sin(q4)+sin(q2+q3).*sin(q5)).^(-1))+J3.*l4.^(-2).* ...
  cos(q5).*(l3.*sin(q2)+l4.*sin(q2+q3)).*sin(q1+(-1).*q4).*(cos(q1) ...
  .*cos(q2+q3).*cos(q4).*cos(q5)+cos(q2+q3).*cos(q5).*sin(q1).*sin( ...
  q4)+sin(q2+q3).*sin(q5)).^(-3).*((dq1+(-1).*dq4).*cos(q1+(-1).*q4) ...
  .*cos(q5).*(l3.*sin(q2)+l4.*sin(q2+q3)).*(cos(q1).*cos(q2+q3).* ...
  cos(q4).*cos(q5)+cos(q2+q3).*cos(q5).*sin(q1).*sin(q4)+sin(q2+q3) ...
  .*sin(q5))+(dq2.*l3.*cos(q2)+(dq2+dq3).*l4.*cos(q2+q3)).*cos(q5).* ...
  sin(q1+(-1).*q4).*(cos(q1).*cos(q2+q3).*cos(q4).*cos(q5)+cos(q2+ ...
  q3).*cos(q5).*sin(q1).*sin(q4)+sin(q2+q3).*sin(q5))+(-1).*dq5.*( ...
  l3.*sin(q2)+l4.*sin(q2+q3)).*sin(q1+(-1).*q4).*sin(q5).*(cos(q1).* ...
  cos(q2+q3).*cos(q4).*cos(q5)+cos(q2+q3).*cos(q5).*sin(q1).*sin(q4) ...
  +sin(q2+q3).*sin(q5))+(-1).*cos(q5).*(l3.*sin(q2)+l4.*sin(q2+q3)) ...
  .*sin(q1+(-1).*q4).*((-1).*dq1.*cos(q2+q3).*cos(q4).*cos(q5).*sin( ...
  q1)+dq4.*cos(q2+q3).*cos(q4).*cos(q5).*sin(q1)+dq5.*cos(q5).*sin( ...
  q2+q3)+(-1).*(dq2+dq3).*cos(q1).*cos(q4).*cos(q5).*sin(q2+q3)+ ...
  dq1.*cos(q1).*cos(q2+q3).*cos(q5).*sin(q4)+(-1).*dq4.*cos(q1).* ...
  cos(q2+q3).*cos(q5).*sin(q4)+(-1).*(dq2+dq3).*cos(q5).*sin(q1).* ...
  sin(q2+q3).*sin(q4)+(dq2+dq3).*cos(q2+q3).*sin(q5)+(-1).*dq5.*cos( ...
  q1).*cos(q2+q3).*cos(q4).*sin(q5)+(-1).*dq5.*cos(q2+q3).*sin(q1).* ...
  sin(q4).*sin(q5)))+l4.^(-1).*lG3.*m3.*cos(q5).*(l3.*sin(q2)+l4.* ...
  sin(q2+q3)).*sin(q1+(-1).*q4).*(cos(q1).*cos(q2+q3).*cos(q4).*cos( ...
  q5)+cos(q2+q3).*cos(q5).*sin(q1).*sin(q4)+sin(q2+q3).*sin(q5)).^( ...
  -1).*((-1).*dq1.*cos(q2+q3).*(l3.*sin(q2)+lG3.*sin(q2+q3))+l4.^( ...
  -1).*lG3.*(cos(q1).*cos(q2+q3).*cos(q4).*cos(q5)+cos(q2+q3).*cos( ...
  q5).*sin(q1).*sin(q4)+sin(q2+q3).*sin(q5)).^(-2).*((dq1+(-1).*dq4) ...
  .*cos(q1+(-1).*q4).*cos(q5).*(l3.*sin(q2)+l4.*sin(q2+q3)).*(cos( ...
  q1).*cos(q2+q3).*cos(q4).*cos(q5)+cos(q2+q3).*cos(q5).*sin(q1).* ...
  sin(q4)+sin(q2+q3).*sin(q5))+(dq2.*l3.*cos(q2)+(dq2+dq3).*l4.*cos( ...
  q2+q3)).*cos(q5).*sin(q1+(-1).*q4).*(cos(q1).*cos(q2+q3).*cos(q4) ...
  .*cos(q5)+cos(q2+q3).*cos(q5).*sin(q1).*sin(q4)+sin(q2+q3).*sin( ...
  q5))+(-1).*dq5.*(l3.*sin(q2)+l4.*sin(q2+q3)).*sin(q1+(-1).*q4).* ...
  sin(q5).*(cos(q1).*cos(q2+q3).*cos(q4).*cos(q5)+cos(q2+q3).*cos( ...
  q5).*sin(q1).*sin(q4)+sin(q2+q3).*sin(q5))+(-1).*cos(q5).*(l3.* ...
  sin(q2)+l4.*sin(q2+q3)).*sin(q1+(-1).*q4).*((-1).*dq1.*cos(q2+q3) ...
  .*cos(q4).*cos(q5).*sin(q1)+dq4.*cos(q2+q3).*cos(q4).*cos(q5).* ...
  sin(q1)+dq5.*cos(q5).*sin(q2+q3)+(-1).*(dq2+dq3).*cos(q1).*cos(q4) ...
  .*cos(q5).*sin(q2+q3)+dq1.*cos(q1).*cos(q2+q3).*cos(q5).*sin(q4)+( ...
  -1).*dq4.*cos(q1).*cos(q2+q3).*cos(q5).*sin(q4)+(-1).*(dq2+dq3).* ...
  cos(q5).*sin(q1).*sin(q2+q3).*sin(q4)+(dq2+dq3).*cos(q2+q3).*sin( ...
  q5)+(-1).*dq5.*cos(q1).*cos(q2+q3).*cos(q4).*sin(q5)+(-1).*dq5.* ...
  cos(q2+q3).*sin(q1).*sin(q4).*sin(q5)))), ...
  dq1.*lG2.^2.*m2.*cos(q2) ...
  .*sin(q2)+l3.*l4.^(-1).*lG3.*m3.*cos(q5).*(l3.*sin(q2)+l4.*sin(q2+ ...
  q3)).*sin(q1+(-1).*q4).*(cos(q1).*cos(q2+q3).*cos(q4).*cos(q5)+ ...
  cos(q2+q3).*cos(q5).*sin(q1).*sin(q4)+sin(q2+q3).*sin(q5)).^(-1).* ...
  ((-1).*dq3.*sin(q3)+(dq2+dq3).*sin(q3)+l4.^(-1).*lG3.*(cos(q1).* ...
  cos(q2).*cos(q4).*cos(q5)+cos(q2).*cos(q5).*sin(q1).*sin(q4)+sin( ...
  q2).*sin(q5)).*(cos(q1).*cos(q2+q3).*cos(q4).*cos(q5)+cos(q2+q3).* ...
  cos(q5).*sin(q1).*sin(q4)+sin(q2+q3).*sin(q5)).^(-2).*((-1).*dq1.* ...
  cos(q2+q3).*cos(q4).*cos(q5).*sin(q1)+dq4.*cos(q2+q3).*cos(q4).* ...
  cos(q5).*sin(q1)+dq5.*cos(q5).*sin(q2+q3)+(-1).*(dq2+dq3).*cos(q1) ...
  .*cos(q4).*cos(q5).*sin(q2+q3)+dq1.*cos(q1).*cos(q2+q3).*cos(q5).* ...
  sin(q4)+(-1).*dq4.*cos(q1).*cos(q2+q3).*cos(q5).*sin(q4)+(-1).*( ...
  dq2+dq3).*cos(q5).*sin(q1).*sin(q2+q3).*sin(q4)+(dq2+dq3).*cos(q2+ ...
  q3).*sin(q5)+(-1).*dq5.*cos(q1).*cos(q2+q3).*cos(q4).*sin(q5)+(-1) ...
  .*dq5.*cos(q2+q3).*sin(q1).*sin(q4).*sin(q5))+l4.^(-1).*lG3.*(cos( ...
  q1).*cos(q2+q3).*cos(q4).*cos(q5)+cos(q2+q3).*cos(q5).*sin(q1).* ...
  sin(q4)+sin(q2+q3).*sin(q5)).^(-1).*((-1).*dq4.*cos(q2).*cos(q4).* ...
  cos(q5).*sin(q1)+(-1).*dq5.*cos(q5).*sin(q2)+dq1.*cos(q2).*cos(q5) ...
  .*sin(q1+(-1).*q4)+dq4.*cos(q1).*cos(q2).*cos(q5).*sin(q4)+dq5.* ...
  cos(q1).*cos(q2).*cos(q4).*sin(q5)+dq5.*cos(q2).*sin(q1).*sin(q4) ...
  .*sin(q5)+dq2.*(cos(q1).*cos(q4).*cos(q5).*sin(q2)+cos(q5).*sin( ...
  q1).*sin(q2).*sin(q4)+(-1).*cos(q2).*sin(q5))))+J3.*l3.*l4.^(-2).* ...
  cos(q5).*(l3.*sin(q2)+l4.*sin(q2+q3)).*sin(q1+(-1).*q4).*(cos(q1) ...
  .*cos(q2+q3).*cos(q4).*cos(q5)+cos(q2+q3).*cos(q5).*sin(q1).*sin( ...
  q4)+sin(q2+q3).*sin(q5)).^(-3).*((cos(q1).*cos(q2).*cos(q4).*cos( ...
  q5)+cos(q2).*cos(q5).*sin(q1).*sin(q4)+sin(q2).*sin(q5)).*((-1).* ...
  dq1.*cos(q2+q3).*cos(q4).*cos(q5).*sin(q1)+dq4.*cos(q2+q3).*cos( ...
  q4).*cos(q5).*sin(q1)+dq5.*cos(q5).*sin(q2+q3)+(-1).*(dq2+dq3).* ...
  cos(q1).*cos(q4).*cos(q5).*sin(q2+q3)+dq1.*cos(q1).*cos(q2+q3).* ...
  cos(q5).*sin(q4)+(-1).*dq4.*cos(q1).*cos(q2+q3).*cos(q5).*sin(q4)+ ...
  (-1).*(dq2+dq3).*cos(q5).*sin(q1).*sin(q2+q3).*sin(q4)+(dq2+dq3).* ...
  cos(q2+q3).*sin(q5)+(-1).*dq5.*cos(q1).*cos(q2+q3).*cos(q4).*sin( ...
  q5)+(-1).*dq5.*cos(q2+q3).*sin(q1).*sin(q4).*sin(q5))+(cos(q1).* ...
  cos(q2+q3).*cos(q4).*cos(q5)+cos(q2+q3).*cos(q5).*sin(q1).*sin(q4) ...
  +sin(q2+q3).*sin(q5)).*((-1).*dq4.*cos(q2).*cos(q4).*cos(q5).*sin( ...
  q1)+(-1).*dq5.*cos(q5).*sin(q2)+dq1.*cos(q2).*cos(q5).*sin(q1+(-1) ...
  .*q4)+dq4.*cos(q1).*cos(q2).*cos(q5).*sin(q4)+dq5.*cos(q1).*cos( ...
  q2).*cos(q4).*sin(q5)+dq5.*cos(q2).*sin(q1).*sin(q4).*sin(q5)+ ...
  dq2.*(cos(q1).*cos(q4).*cos(q5).*sin(q2)+cos(q5).*sin(q1).*sin(q2) ...
  .*sin(q4)+(-1).*cos(q2).*sin(q5))))+dq1.*l3.*m3.*(l3.*sin(q2)+ ...
  lG3.*sin(q2+q3)).*(sin(q3).*sin(q2+q3)+cos(q2+q3).*(cos(q3)+(-1).* ...
  l4.^(-1).*lG3.*(cos(q1).*cos(q2).*cos(q4).*cos(q5)+cos(q2).*cos( ...
  q5).*sin(q1).*sin(q4)+sin(q2).*sin(q5)).*(cos(q1).*cos(q2+q3).* ...
  cos(q4).*cos(q5)+cos(q2+q3).*cos(q5).*sin(q1).*sin(q4)+sin(q2+q3) ...
  .*sin(q5)).^(-1)));
  (-1/2).*dq1.*lG2.^2.*m2.*sin(2.*q2)+l3.*m3.* ...
  sin(q3).*((-1).*dq1.*sin(q2+q3).*(l3.*sin(q2)+lG3.*sin(q2+q3))+( ...
  -1).*(dq2+dq3).*l4.^(-1).*lG3.*cos(q5).*(l3.*sin(q2)+l4.*sin(q2+ ...
  q3)).*sin(q1+(-1).*q4).*(cos(q1).*cos(q2+q3).*cos(q4).*cos(q5)+ ...
  cos(q2+q3).*cos(q5).*sin(q1).*sin(q4)+sin(q2+q3).*sin(q5)).^(-1))+ ...
  (-1).*J3.*l3.*l4.^(-2).*(cos(q1).*cos(q2).*cos(q4).*cos(q5)+cos( ...
  q2).*cos(q5).*sin(q1).*sin(q4)+sin(q2).*sin(q5)).*(cos(q1).*cos( ...
  q2+q3).*cos(q4).*cos(q5)+cos(q2+q3).*cos(q5).*sin(q1).*sin(q4)+ ...
  sin(q2+q3).*sin(q5)).^(-3).*((dq1+(-1).*dq4).*cos(q1+(-1).*q4).* ...
  cos(q5).*(l3.*sin(q2)+l4.*sin(q2+q3)).*(cos(q1).*cos(q2+q3).*cos( ...
  q4).*cos(q5)+cos(q2+q3).*cos(q5).*sin(q1).*sin(q4)+sin(q2+q3).* ...
  sin(q5))+(dq2.*l3.*cos(q2)+(dq2+dq3).*l4.*cos(q2+q3)).*cos(q5).* ...
  sin(q1+(-1).*q4).*(cos(q1).*cos(q2+q3).*cos(q4).*cos(q5)+cos(q2+ ...
  q3).*cos(q5).*sin(q1).*sin(q4)+sin(q2+q3).*sin(q5))+(-1).*dq5.*( ...
  l3.*sin(q2)+l4.*sin(q2+q3)).*sin(q1+(-1).*q4).*sin(q5).*(cos(q1).* ...
  cos(q2+q3).*cos(q4).*cos(q5)+cos(q2+q3).*cos(q5).*sin(q1).*sin(q4) ...
  +sin(q2+q3).*sin(q5))+(-1).*cos(q5).*(l3.*sin(q2)+l4.*sin(q2+q3)) ...
  .*sin(q1+(-1).*q4).*((-1).*dq1.*cos(q2+q3).*cos(q4).*cos(q5).*sin( ...
  q1)+dq4.*cos(q2+q3).*cos(q4).*cos(q5).*sin(q1)+dq5.*cos(q5).*sin( ...
  q2+q3)+(-1).*(dq2+dq3).*cos(q1).*cos(q4).*cos(q5).*sin(q2+q3)+ ...
  dq1.*cos(q1).*cos(q2+q3).*cos(q5).*sin(q4)+(-1).*dq4.*cos(q1).* ...
  cos(q2+q3).*cos(q5).*sin(q4)+(-1).*(dq2+dq3).*cos(q5).*sin(q1).* ...
  sin(q2+q3).*sin(q4)+(dq2+dq3).*cos(q2+q3).*sin(q5)+(-1).*dq5.*cos( ...
  q1).*cos(q2+q3).*cos(q4).*sin(q5)+(-1).*dq5.*cos(q2+q3).*sin(q1).* ...
  sin(q4).*sin(q5)))+l3.*m3.*(cos(q3)+(-1).*l4.^(-1).*lG3.*(cos(q1) ...
  .*cos(q2).*cos(q4).*cos(q5)+cos(q2).*cos(q5).*sin(q1).*sin(q4)+ ...
  sin(q2).*sin(q5)).*(cos(q1).*cos(q2+q3).*cos(q4).*cos(q5)+cos(q2+ ...
  q3).*cos(q5).*sin(q1).*sin(q4)+sin(q2+q3).*sin(q5)).^(-1)).*((-1) ...
  .*dq1.*cos(q2+q3).*(l3.*sin(q2)+lG3.*sin(q2+q3))+l4.^(-1).*lG3.*( ...
  cos(q1).*cos(q2+q3).*cos(q4).*cos(q5)+cos(q2+q3).*cos(q5).*sin(q1) ...
  .*sin(q4)+sin(q2+q3).*sin(q5)).^(-2).*((dq1+(-1).*dq4).*cos(q1+( ...
  -1).*q4).*cos(q5).*(l3.*sin(q2)+l4.*sin(q2+q3)).*(cos(q1).*cos(q2+ ...
  q3).*cos(q4).*cos(q5)+cos(q2+q3).*cos(q5).*sin(q1).*sin(q4)+sin( ...
  q2+q3).*sin(q5))+(dq2.*l3.*cos(q2)+(dq2+dq3).*l4.*cos(q2+q3)).* ...
  cos(q5).*sin(q1+(-1).*q4).*(cos(q1).*cos(q2+q3).*cos(q4).*cos(q5)+ ...
  cos(q2+q3).*cos(q5).*sin(q1).*sin(q4)+sin(q2+q3).*sin(q5))+(-1).* ...
  dq5.*(l3.*sin(q2)+l4.*sin(q2+q3)).*sin(q1+(-1).*q4).*sin(q5).*( ...
  cos(q1).*cos(q2+q3).*cos(q4).*cos(q5)+cos(q2+q3).*cos(q5).*sin(q1) ...
  .*sin(q4)+sin(q2+q3).*sin(q5))+(-1).*cos(q5).*(l3.*sin(q2)+l4.* ...
  sin(q2+q3)).*sin(q1+(-1).*q4).*((-1).*dq1.*cos(q2+q3).*cos(q4).* ...
  cos(q5).*sin(q1)+dq4.*cos(q2+q3).*cos(q4).*cos(q5).*sin(q1)+dq5.* ...
  cos(q5).*sin(q2+q3)+(-1).*(dq2+dq3).*cos(q1).*cos(q4).*cos(q5).* ...
  sin(q2+q3)+dq1.*cos(q1).*cos(q2+q3).*cos(q5).*sin(q4)+(-1).*dq4.* ...
  cos(q1).*cos(q2+q3).*cos(q5).*sin(q4)+(-1).*(dq2+dq3).*cos(q5).* ...
  sin(q1).*sin(q2+q3).*sin(q4)+(dq2+dq3).*cos(q2+q3).*sin(q5)+(-1).* ...
  dq5.*cos(q1).*cos(q2+q3).*cos(q4).*sin(q5)+(-1).*dq5.*cos(q2+q3).* ...
  sin(q1).*sin(q4).*sin(q5)))), ...
  l3.^2.*(m3.*(cos(q3)+(-1).*l4.^(-1).* ...
  lG3.*(cos(q1).*cos(q2).*cos(q4).*cos(q5)+cos(q2).*cos(q5).*sin(q1) ...
  .*sin(q4)+sin(q2).*sin(q5)).*(cos(q1).*cos(q2+q3).*cos(q4).*cos( ...
  q5)+cos(q2+q3).*cos(q5).*sin(q1).*sin(q4)+sin(q2+q3).*sin(q5)).^( ...
  -1)).*((-1).*dq3.*sin(q3)+(dq2+dq3).*sin(q3)+l4.^(-1).*lG3.*(cos( ...
  q1).*cos(q2).*cos(q4).*cos(q5)+cos(q2).*cos(q5).*sin(q1).*sin(q4)+ ...
  sin(q2).*sin(q5)).*(cos(q1).*cos(q2+q3).*cos(q4).*cos(q5)+cos(q2+ ...
  q3).*cos(q5).*sin(q1).*sin(q4)+sin(q2+q3).*sin(q5)).^(-2).*((-1).* ...
  dq1.*cos(q2+q3).*cos(q4).*cos(q5).*sin(q1)+dq4.*cos(q2+q3).*cos( ...
  q4).*cos(q5).*sin(q1)+dq5.*cos(q5).*sin(q2+q3)+(-1).*(dq2+dq3).* ...
  cos(q1).*cos(q4).*cos(q5).*sin(q2+q3)+dq1.*cos(q1).*cos(q2+q3).* ...
  cos(q5).*sin(q4)+(-1).*dq4.*cos(q1).*cos(q2+q3).*cos(q5).*sin(q4)+ ...
  (-1).*(dq2+dq3).*cos(q5).*sin(q1).*sin(q2+q3).*sin(q4)+(dq2+dq3).* ...
  cos(q2+q3).*sin(q5)+(-1).*dq5.*cos(q1).*cos(q2+q3).*cos(q4).*sin( ...
  q5)+(-1).*dq5.*cos(q2+q3).*sin(q1).*sin(q4).*sin(q5))+l4.^(-1).* ...
  lG3.*(cos(q1).*cos(q2+q3).*cos(q4).*cos(q5)+cos(q2+q3).*cos(q5).* ...
  sin(q1).*sin(q4)+sin(q2+q3).*sin(q5)).^(-1).*((-1).*dq4.*cos(q2).* ...
  cos(q4).*cos(q5).*sin(q1)+(-1).*dq5.*cos(q5).*sin(q2)+dq1.*cos(q2) ...
  .*cos(q5).*sin(q1+(-1).*q4)+dq4.*cos(q1).*cos(q2).*cos(q5).*sin( ...
  q4)+dq5.*cos(q1).*cos(q2).*cos(q4).*sin(q5)+dq5.*cos(q2).*sin(q1) ...
  .*sin(q4).*sin(q5)+dq2.*(cos(q1).*cos(q4).*cos(q5).*sin(q2)+cos( ...
  q5).*sin(q1).*sin(q2).*sin(q4)+(-1).*cos(q2).*sin(q5))))+(-1).* ...
  J3.*l4.^(-2).*(cos(q1).*cos(q2).*cos(q4).*cos(q5)+cos(q2).*cos(q5) ...
  .*sin(q1).*sin(q4)+sin(q2).*sin(q5)).*(cos(q1).*cos(q2+q3).*cos( ...
  q4).*cos(q5)+cos(q2+q3).*cos(q5).*sin(q1).*sin(q4)+sin(q2+q3).* ...
  sin(q5)).^(-3).*((cos(q1).*cos(q2).*cos(q4).*cos(q5)+cos(q2).*cos( ...
  q5).*sin(q1).*sin(q4)+sin(q2).*sin(q5)).*((-1).*dq1.*cos(q2+q3).* ...
  cos(q4).*cos(q5).*sin(q1)+dq4.*cos(q2+q3).*cos(q4).*cos(q5).*sin( ...
  q1)+dq5.*cos(q5).*sin(q2+q3)+(-1).*(dq2+dq3).*cos(q1).*cos(q4).* ...
  cos(q5).*sin(q2+q3)+dq1.*cos(q1).*cos(q2+q3).*cos(q5).*sin(q4)+( ...
  -1).*dq4.*cos(q1).*cos(q2+q3).*cos(q5).*sin(q4)+(-1).*(dq2+dq3).* ...
  cos(q5).*sin(q1).*sin(q2+q3).*sin(q4)+(dq2+dq3).*cos(q2+q3).*sin( ...
  q5)+(-1).*dq5.*cos(q1).*cos(q2+q3).*cos(q4).*sin(q5)+(-1).*dq5.* ...
  cos(q2+q3).*sin(q1).*sin(q4).*sin(q5))+(cos(q1).*cos(q2+q3).*cos( ...
  q4).*cos(q5)+cos(q2+q3).*cos(q5).*sin(q1).*sin(q4)+sin(q2+q3).* ...
  sin(q5)).*((-1).*dq4.*cos(q2).*cos(q4).*cos(q5).*sin(q1)+(-1).* ...
  dq5.*cos(q5).*sin(q2)+dq1.*cos(q2).*cos(q5).*sin(q1+(-1).*q4)+ ...
  dq4.*cos(q1).*cos(q2).*cos(q5).*sin(q4)+dq5.*cos(q1).*cos(q2).* ...
  cos(q4).*sin(q5)+dq5.*cos(q2).*sin(q1).*sin(q4).*sin(q5)+dq2.*( ...
  cos(q1).*cos(q4).*cos(q5).*sin(q2)+cos(q5).*sin(q1).*sin(q2).*sin( ...
  q4)+(-1).*cos(q2).*sin(q5))))+m3.*sin(q3).*(dq3.*cos(q3)+((-1).* ...
  dq2+(-1).*dq3).*(cos(q3)+(-1).*l4.^(-1).*lG3.*(cos(q1).*cos(q2).* ...
  cos(q4).*cos(q5)+cos(q2).*cos(q5).*sin(q1).*sin(q4)+sin(q2).*sin( ...
  q5)).*(cos(q1).*cos(q2+q3).*cos(q4).*cos(q5)+cos(q2+q3).*cos(q5).* ...
  sin(q1).*sin(q4)+sin(q2+q3).*sin(q5)).^(-1))))];
