function memberFEF = computedFEF(w,L)
%computes the fixed end forces for each element
        memberFEF = [-w(1,1)*L/2, -w(1,2)*L/2, 0,0, 0, -w(1,2)*L^2/12,...
        -w(1,1)*L/2, -w(1,2)*L/2, 0,0, 0, w(1,2)*L^2/12];
end 





