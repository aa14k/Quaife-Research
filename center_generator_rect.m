% This functions generates random centers that are picked such that
% No shape intersects another shape

function [X,Y] = center_generator_rect(r,N,xl,xr,yb,yt)

X = xl + (xr-xl)*rand(N,1);
Y = yb + (yt-yb)*rand(N,1); 
X(1) = xl + (xr-xl)*rand();
Y(1) = yb + (yt-yb)*rand();
i = 2;
while(i<(N+1))
    x = xl + (xr-xl)*rand();
    y = yb + (yt-yb)*rand();
    if 1
        distance = zeros(i-1,1);
        for j=1:i
            distance(j) = sqrt((x-X(j))^2 + (y-Y(j))^2);
        end
        disTest = distance < r*2;
        if sum(disTest) == 0
            X(i) = x;
            Y(i) = y;
            i = i + 1;
        end
    end
end

end
