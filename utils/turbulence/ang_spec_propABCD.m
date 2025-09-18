function [x2 y2 Uout] ...
    = ang_spec_propABCD(Uin, wvl, d1, d2, ABCD)
% function [x2 y2 Uout] ...
%     = ang_spec_propABCD(Uin, wwl, d1, d2, ABCD)

    N = size(Uin,1);   % assume square grid
    k = 2*pi/wvl;    % optical wavevector
    % source-plane coordinates
    [x1 y1] = meshgrid((-N/2 : 1 : N/2 - 1) * d1);
    r1sq = x1.^2 + y1.^2;
    % spatial frequencies (of source plane)
    df1 = 1 / (N*d1);
    [fX fY] = meshgrid((-N/2 : 1 : N/2 - 1) * df1);
    fsq = fX.^2 + fY.^2;
    % scaling parameter
    m = d2/d1;
    % observation-plane coordinates
    [x2 y2] = meshgrid((-N/2 : 1 : N/2 - 1) * d2);
    r2sq = x2.^2 + y2.^2;
    % optical system matrix
    A = ABCD(1,1); B = ABCD(1,2); D = ABCD(2,2);
    % quadratic phase factors
    Q1 = exp(i*pi/(wvl*B)*(A-m)*r1sq);
    Q2 = exp(-i*pi*wvl*B/m*fsq);
    Q3 = exp(i*pi/(wvl*B)*(D-1/m)*r2sq);
    % compute the propagated field
    Uout = Q3.* ift2(Q2 .* ft2(Q1 .* Uin / m, d1), df1);