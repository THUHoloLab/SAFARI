function U = genAiry(X,Y,w0,x0,y0,a)
% =========================================================================
% Generate Airy beams.
% -------------------------------------------------------------------------
% Input:    - X / Y   : 2D coordinates (mm).
%           - w0      : Scaling factor.
%           - x0 / y0 : Center location (mm).
%           - a       : Exponential truncation factor.
% Output:   - U       : Complex amplitude of the wavefront.
% =========================================================================
U = airy((x0 - X)/w0) .* airy((y0 - Y)/w0) ...
    .* exp(a*(x0 - X)/w0) .* exp(a*(y0 - Y)/w0);

end