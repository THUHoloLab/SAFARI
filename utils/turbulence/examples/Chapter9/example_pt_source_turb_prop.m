% example_pt_source_turb_prop.m

l0 = 0;     % inner scale [m]
L0 = inf;   % outer scale [m]

zt = [0 z];  % propagation plane locations
Delta_z = zt(2:n) - zt(1:n-1);    % propagation distances
% grid spacings
alpha = zt / zt(n);
delta = (1-alpha) * delta1 + alpha * deltan;

% initialize array for phase screens
phz = zeros(N, N, n);
nreals = 20;    % number of random realizations
% initialize arrays for propagated fields,
% aperture mask, and MCF
Uout = zeros(N);
mask = circ(xn/D2, yn/D2, 1);
MCF2 = zeros(N);
sg = repmat(sg, [1 1 n]);
for idxreal = 1 : nreals     % loop over realizations
    idxreal
    % loop over screens
    for idxscr = 1 : 1 : n
        [phz_lo phz_hi] ...
            = ft_sh_phase_screen ...
            (r0scrn(idxscr), N, delta(idxscr), L0, l0);
        phz(:,:,idxscr) = phz_lo + phz_hi;
    end
    % simulate turbulent propagation
    [xn yn Uout] = ang_spec_multi_prop(pt, wvl, ....
        delta1, deltan, z, sg.*exp(i*phz));
    % collimate the beam
    Uout = Uout .* exp(-i*pi/(wvl*R)*(xn.^2+yn.^2));
    % accumulate realizations of the MCF
    MCF2 = MCF2 + corr2_ft(Uout, Uout, mask, deltan);
end
% modulus of the complex degree of coherence
MCDOC2 = abs(MCF2) / (MCF2(N/2+1,N/2+1));