function imgfit = legendrepolyfit(xx, yy, orders, coefs)

o1 = orders(1);
o2 = orders(2);

imgfit = zeros(size(xx));
for i = 0:o1
    for j = 0:o2
        
        index = i*(o2+1)+j+1;
        imgfit = imgfit + coefs(index) * myLegendreP(i,xx).*myLegendreP(j,yy);
    end
end

end

