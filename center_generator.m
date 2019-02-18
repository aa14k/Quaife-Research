% This functions generates random centers that are picked such that
% No shape intersects another shape

function [X,Y] = center_generator(r,N,xl,xr,yb,yt)

X = zeros(N,1);
Y = zeros(N,1); 
X(1) = xl + (xr-xl)*rand();
Y(1) = yb + (yt-yb)*rand();
i = 2;
l = 1;
while(i<(N+1))
    x = xl + (xr-xl)*rand();
    y = yb + (yt-yb)*rand();
    distance = zeros(i-1,1);
    for j=1:i
        distance(j) = sqrt((x-X(j))^2 + (y-Y(j))^2);
    end
    disTest = distance < r;
    if sum(disTest) == 0
        X(i) = x;
        Y(i) = y;
        i = i + 1;
    end
    l = l + 1;
end

end

