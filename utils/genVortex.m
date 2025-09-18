function phase = genVortex(x,y,c,order)

phase = atan((y-c(2))./(x-c(1)));

mask  = (x - c(1) < 0) & (y - c(2) >= 0);
index = find(mask > 0);
phase(index) = phase(index) + pi;

mask = (x - c(1) < 0) & (y - c(2) < 0);
index = find(mask > 0);
phase(index) = phase(index) - pi;

phase = phase * order;

end

