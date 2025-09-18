function y = myLaguerre(p,l,x)
% =========================================================================
% Calculate the Laguerre functions.
% -------------------------------------------------------------------------
% Input:    - x : Input.
%           - l : Azimuthal index.
%           - p : Radial index.
% Output:   - y : Function value.
% =========================================================================

y = zeros(p+1,1);
if p == 0
    y = 1;
else
for m = 0:p
    y(p+1-m) = ((-1).^m.*(factorial(p+l)))./(factorial(p-m).*factorial(l+m).*factorial(m));
end
end
y = polyval(y,x);