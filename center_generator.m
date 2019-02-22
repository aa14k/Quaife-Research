% This functions generates random centers that are picked such that
% No shape intersects another shape

function [X,Y] = center_generator(r,N,xl,xr,yb,yt,R)
R = R-r;
X = xl + (xr-xl)*rand(N,1);
Y = yb + (yt-yb)*rand(N,1); 
X(1) = xl + (xr-xl)*rand();
Y(1) = yb + (yt-yb)*rand();
while X(1)^2 + Y(1)^2 > R^2
    X(1) = xl + (xr-xl)*rand();
    Y(1) = yb + (yt-yb)*rand();
end
i = 2;
while(i<(N+1))
    x = xl + (xr-xl)*rand();
    y = yb + (yt-yb)*rand();
    if x^2 + y^2 < R^2
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
    end
end

end

