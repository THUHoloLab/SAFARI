function c = sinebow(n)
h = linspace(0,1,n);
h = h + 1/2;
h = h * (-1);
r = sin(pi*h);
g = sin(pi*(h+1/3));
b = sin(pi*(h+2/3));
c = [r;g;b]';
c = c.^2;
end