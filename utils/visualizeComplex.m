function [img,cbarimg] = visualizeComplex(cimg, cmap, amp_range, mode, reverse)

amp = abs(cimg);
pha = angle(cimg);

amin = amp_range(1);
amax = amp_range(2);

amp_norm = (amp-amin)/(amax-amin);

b = 1;

ncmap = length(cmap);
img = zeros(size(cimg,1),size(cimg,2),3);

for i = 1:size(cimg,1)
    for j = 1:size(cimg,2)
        a = amp_norm(i,j);
        if reverse
            a = 1-a;
        end
        if strcmpi(mode,'hsv')
            img(i,j,:) = cmap(1+round((ncmap-1)/2/pi*(pha(i,j)+pi)),:) .* a;
        elseif strcmpi(mode,'hsl')
            if a > 1/2
                w = (2*(a-1/2))^b;
                img(i,j,:) = cmap(1+round((ncmap-1)/2/pi*(pha(i,j)+pi)),:) .* (1-w) + [1,1,1] .* w;
            else
                w = (2*(1/2-a))^b;
                img(i,j,:) = cmap(1+round((ncmap-1)/2/pi*(pha(i,j)+pi)),:) .* (1-w) + [0,0,0] .* w;
            end
        end
        
    end
end

n = 256;
cbarimg = nan(n,n,3);
x = linspace(-1,1,n);
y = x;
[X,Y] = meshgrid(x,y);
[theta,rho] = cart2pol(X,Y);
for i = 1:n
    for j = 1:n
        if rho(i,j) <= 1
            a = rho(i,j);
            if reverse
                a = 1-a;
            end
            if strcmpi(mode,'hsv')
                cbarimg(i,j,:) = cmap(1+round((ncmap-1)/2/pi*(theta(i,j)+pi)),:) .* a;
            elseif strcmpi(mode,'hsl')
                if a > 1/2
                    w = (2*(a-1/2))^b;
                    cbarimg(i,j,:) = cmap(1+round((ncmap-1)/2/pi*(theta(i,j)+pi)),:) .* (1-w) + [1,1,1] .* w;
                else
                    w = (2*(1/2-a))^b;
                    cbarimg(i,j,:) = cmap(1+round((ncmap-1)/2/pi*(theta(i,j)+pi)),:) .* (1-w) + [0,0,0] .* w;
                end
            end
        end
    end
end

end