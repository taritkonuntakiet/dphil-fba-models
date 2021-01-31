function [dBCdya, dBCdyb] = flowering4bcjac (ya, yb, v, ph)
% flowerin4bcjac provides jacobian of BCs.

dBCdya = eye(15);       % ya(i) only appears in the i-th eqn. 
dBCdyb = - eye(15);     % -yb(i) only appears in the i-th eqn.

