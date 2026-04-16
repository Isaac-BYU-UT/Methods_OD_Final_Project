function [uCross] = crossProductEquivalent(u)
% crossProductEquivalent : Outputs the cross-product-equivalent matrix uCross 
%                           such that for arbitrary 3-by-1 vectors u and v,
%                           cross(u,v) = uCross*v.
%
%
% INPUTS
%
% u---------- 3-by-1 vector
%
%
% OUTPUTS
%
% uCross----- 3-by-3 skew-symmetric cross-product equivalent matrix
%
%+------------------------------------------------------------------------------+
% References: Dr. Humphreys' main.pdf (pg. 7)
% (https://gitlab.com/todd.humphreys/ar-book/-/blob/master/main.pdf?ref_type=heads)
%
%
% Author: Isaac Sorensen (Validated on 1/23/2026)
%+==============================================================================+

% Let's gaurantee that the input is a column vector:

u = u(:);

u1 = u(1);
u2 = u(2);
u3 = u(3);
uCross = [  0  -u3 u2;...
            u3  0 -u1;...
            -u2 u1 0];

end