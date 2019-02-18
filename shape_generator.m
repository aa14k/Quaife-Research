%hyper parameters
th = 0:pi/50:2*pi;
w = 0.1; % closer to 0 the more "circular" the shape will be
r0 = 1; %radius
x0 = 3; % center of blob for x axis
y0 = 2; % center of blob for y axis

%random weights intialized
t1 = 2*pi*.468;
t2 = 2*pi*0.9;
t3 = 2*pi*.464;

w1 = 0.225;
w2 = 0.159;
w3 = 0.225;

% our "blob"
r=r0-w*(w1*sin(th+t1)+w2*sin(2*th+t2)+w3*sin(3*th+t3)); %creates a "blob"
r1 = r0 + r0*w;
% plot in cartesian
[x,y] = pol2cart(th,r);
[x1,y1] = pol2cart(th,r1);
x = x + x0;
y = y + y0;
x1 = x1 + x0;
y1 = y1 + y0;
plot(x,y)
hold on
plot(x1,y1)
scatter(x0,y0)
hold off
