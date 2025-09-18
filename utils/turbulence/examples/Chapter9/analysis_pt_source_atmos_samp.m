% analysis_pt_source_atmos_samp.m

c = 2;
D1p = D1 + c*wvl*Dz/r0sw;
D2p = D2 + c*wvl*Dz/r0sw;

delta1 = linspace(0, 1.1*wvl*Dz/D2p, 100);
deltan = linspace(0, 1.1*wvl*Dz/D1p, 100);
% constraint 1
deltan_max = -D2p/D1p*delta1 + wvl*Dz/D1p;
% constraint 3
d2min3 = (1+Dz/R)*delta1 - wvl*Dz/D1p;
d2max3 = (1+Dz/R)*delta1 + wvl*Dz/D1p;
[delta1 deltan] = meshgrid(delta1, deltan);
% constraint 2
N2 = (wvl * Dz + D1p*deltan + D2p*delta1) ...
    ./ (2 * delta1 .* deltan);

% constraint 4
d1 = 10e-3;
d2 = 10e-3;
N = 512;
d1*d2 * N / wvl
zmax = min([d1 d2])^2 * N / wvl
nmin = ceil(Dz / zmax) + 1