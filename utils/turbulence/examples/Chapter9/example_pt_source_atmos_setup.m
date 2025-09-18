% example_pt_source_atmos_setup.m

% determine geometry
D2 = 0.5;   % diameter of the observation aperture [m]
wvl = 1e-6;  % optical wavelength [m]
k = 2*pi / wvl; % optical wavenumber [rad/m]
Dz = 50e3;     % propagation distance [m]

% use sinc to model pt source
DROI = 4 * D2;  % diam of obs-plane region of interest [m]
D1 = wvl*Dz / DROI;    % width of central lobe [m]
R = Dz; % wavefront radius of curvature [m]

% atmospheric properties
Cn2 = 1e-16;    % structure parameter [m^-2/3], constant
% SW and PW coherence diameters [m]
r0sw = (0.423 * k^2 * Cn2 * 3/8 * Dz)^(-3/5);
r0pw = (0.423 * k^2 * Cn2 * Dz)^(-3/5);
p = linspace(0, Dz, 1e3);
% log-amplitude variance
rytov = 0.563 * k^(7/6) * sum(Cn2 * (1-p/Dz).^(5/6) ...
    .* p.^(5/6) * (p(2)-p(1)));

% screen properties
nscr = 11; % number of screens
A = zeros(2, nscr); % matrix
alpha = (0:nscr-1) / (nscr-1);
A(1,:) = alpha.^(5/3);
A(2,:) = (1 - alpha).^(5/6) .* alpha.^(5/6);
b = [r0sw.^(-5/3); rytov/1.33*(k/Dz)^(5/6)];
% initial guess
x0 = (nscr/3*r0sw * ones(nscr, 1)).^(-5/3);
% objective function
fun = @(X) sum((A*X(:) - b).^2);
% constraints
x1 = zeros(nscr, 1);
rmax = 0.1; % maximum Rytov number per partial prop
x2 = rmax/1.33*(k/Dz)^(5/6) ./ A(2,:);
x2(A(2,:)==0) = 50^(-5/3)
[X,fval,exitflag,output] ...
    = fmincon(fun,x0,[],[],[],[],x1,x2)
% check screen r0s
r0scrn = X.^(-3/5)
r0scrn(isinf(r0scrn)) = 1e6;
% check resulting r0sw & rytov
bp = A*X(:); [bp(1)^(-3/5) bp(2)*1.33*(Dz/k)^(5/6)]
[r0sw rytov]