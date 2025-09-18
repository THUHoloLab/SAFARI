function [A, xx, yy, rects] = legendrebasis(img, orders, n_rect)

n1 = size(img,1);
n2 = size(img,2);

o1 = orders(1);
o2 = orders(2);

rects = nan(n_rect,4);

for i = 1:n_rect
    fig = figure;
    [~,rect_tmp] = imcrop((img - min(img(:)))/(max(img(:)) - min(img(:))));
    close(fig);
    rects(i,:) = rect_tmp;
end

c1 = linspace(0,1,n1);
c2 = linspace(0,1,n2);

[xx,yy] = meshgrid(c2,c1);

xxx = [];
yyy = [];

for i = 1:n_rect
    xx_crop = imcrop(xx, rects(i,:));
    yy_crop = imcrop(yy, rects(i,:));
    xxx = [xxx; xx_crop(:)];
    yyy = [yyy; yy_crop(:)];
end



n_coef = (o1+1)*(o2+1);
A = nan(length(xxx),n_coef);

for i = 0:o1
    for j = 0:o2
        
        index = i*(o2+1)+j+1;
%         A(:,index) = legendreP(i,xxx).*legendreP(j,yyy);
        A(:,index) = myLegendreP(i,xxx).*myLegendreP(j,yyy);
    end
end

end

