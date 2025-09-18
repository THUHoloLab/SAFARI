function [out] = DMD_Hol(U, X, Y, gx, gy, rot, weight)
    if rot == 0
        A = abs(U); % correction factors
        A = (A*weight)./max(A(:));
        A = asin(A)/pi; % amplitude term for hologram
        Phi = angle(U); Phi = Phi/pi; % phase term for hologram
        h = 0.5+0.5*sign(cos(2*pi*(gx*X+gy*Y)+pi*Phi)-cos(pi*A)); % generates complex amplitude hologram
        % h = 0.5+0.5*sign(cos(2*pi)-cos(pi*A)); % generates amplitude only hologram
        % h = 0.5+0.5*sign(cos(2*pi*(gx*X+gy*Y)+pi*Phi)); % generates phase only hologram

        h(h==0.5) = 0; % binarises hologram
        normh = h/max(h(:));
        normh = normh';
    else
        U0 = U;
        U = imrotate(U,rot);
        U = U(size(U,1)/2 - size(U0,1)/2:size(U,1)/2 + size(U0,1)/2 - 1, ...
              size(U,2)/2 - size(U0,2)/2:size(U,2)/2 + size(U0,2)/2 - 1);
        A = abs(U); % correction factors
        A = (A*weight)./max(A(:));
        A = asin(A)/pi; % amplitude term for hologram
        Phi = angle(U); Phi = Phi/pi; % phase term for hologram
        h = 0.5+0.5*xsign(cos(2*pi*(gx*X+gy*Y)+pi*Phi)-cos(pi*A)); % generates hologram
        % h = 0.5+0.5*sign(cos(2*pi)-cos(pi*A)); % generates amplitude only hologram
        % h = 0.5+0.5*sign(cos(2*pi*(gx*X+gy*Y)+pi*Phi)); % generates phase only hologram

        h(h==0.5) = 0; % binarises hologram
        normh = h/max(h(:));
        normh = normh';
    end
    out = normh;

end