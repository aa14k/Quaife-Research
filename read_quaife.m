% Read velocity file
fileName = 'circles17EulerVelocities.bin';
[nx,ny,eulerX,eulerY,u,v] = ...
    loadEulerVelocities(fileName);

% Read trap file
Xcenters_file = matfile('Xcenters.mat');
Ycenters_file = matfile('Ycenters.mat');
radii_file = matfile('radii.mat');

% Load data into arrays
Xcenters = Xcenters_file.Xcenters;
Ycenters = Ycenters_file.Ycenters;

radii = radii_file.radii;
centers = Xcenters + 1i*Ycenters;
max_rad = max(radii);

% Discretization Interval
nx = 500; %number of points in the x direction
ny = 500; %number of points in the y direction

% Set bounding shape
xmin = 1;
xmax = 2;
ymin = 1;
ymax = 2;

% Plot only a subsection of the circles
sL = real(centers) > xmin - 0.3;
sR = real(centers) < xmax + 0.3;
sU = imag(centers) > ymin - 0.3;
sD = imag(centers) < ymax + 0.3;

centers = centers( sL & sR & sU & sD );
radii = radii( sL & sR & sU & sD );

N = numel(centers);

% Create interpolating objects that can be queried at given points
Fx = griddedInterpolant( eulerX', eulerY', u', 'linear' );
Fy = griddedInterpolant( eulerX', eulerY', v', 'linear' );

% Definte initial walker positions
[x,y] = meshgrid(linspace(xmin,xmax,nx),linspace(ymin,ymax,ny));
init_pts = x(:) + 1i*y(:);% initial location of the random walks
npts = numel(init_pts);% number of initial points

nruns = 100;% number of runs per initial location
mfpt = zeros(npts,1);% initialize our mean first passage time vector

p1 = 1;
p2 = 0;

dt = 0.1;% our step size for the random walker
for k = 1:numel(init_pts)%iterate over each of the nx by ny inital pts
    disp(init_pts(k))
      
    % Initialize walkers
    z = init_pts(k)*ones(nruns,1); 
    zold = z;
    znew = zold;
     
    inside = zeros(size(z));
    for i=1:nruns
        inside(i) = any( abs( znew(i) - centers ) < radii );
    end      


    fpt = zeros(nruns,1);
    % min_dist.m
    while any(~inside) 
        s = find(~inside);

        fpt(s) = fpt(s) + dt;
        omega = 2*pi*rand(nruns,1);
        alpha = rand(nruns, 1);

        xinterp = Fx( real(zold(s)), imag(zold(s)) );
        yinterp = Fy( real(zold(s)), imag(zold(s)) );

        znew(s) = zold(s) + p1*dt*exp(1i*omega(s)) .* sqrt(-dt*2*log(alpha(s))) + ...
                            p2*dt*(xinterp + 1i*yinterp);

                        
        inside = zeros(size(z));
        for i=1:nruns
            inside(i) = any( abs( znew(i) - centers ) < radii );
        end       

        % Check top bound
        if any((imag(znew) > ymax))
            s = find((imag(znew) > ymax));    
            %znew(s) = real(znew(s)) + 1i*ymax; % Reflecting
            inside(s) = 1; % Absorbing
        end

        % Check bottom bound
        if any((imag(znew) < ymin))
            s = find((imag(znew) < ymin));    
            %znew(s) = real(znew(s)) + 1i*ymin; % Reflecting
            inside(s) = 1; % Absorbing 
        end

        % Check left bound
        if any((real(znew) < xmin))
            s = find((real(znew) < xmin));    
            %znew(s) = xmin + 1i*imag(znew(s)); % Reflecting
            inside(s) = 1; % Absorbing
        end

        % Check right bound
        if any((real(znew) > xmax))
            s = find((real(znew) > xmax));    
            %znew(s) = xmax + 1i*imag(znew(s)); % Reflecting
            inside(s) = 1; % Absorbing
        end

        zold = znew;
    end

    mfpt(k) = mean(fpt); % stores the mean first passage of the random walker
end

% plot circles
figure(1);
clf;
mfpt = reshape(mfpt,ny,nx);
ss = find(mfpt < 0.01*dt);
mfpt(ss) = 1/0;

surf(x,y,mfpt);
colorbar
view(2); shading interp;

axis equal; hold on

for i=1:N
    th = 0:pi/50:2*pi;
    xunit = radii(i) * cos(th) + real(centers(i));
    yunit = radii(i) * sin(th) + imag(centers(i));
    h = fill(xunit, yunit, 'k');
    set(h,'edgecolor',[.5 .5 .5])
    set(h,'facecolor',[.75 .75 .75])
end


xlim([xmin xmax])
ylim([ymin ymax])
hold off

function [nx,ny,eulerX,eulerY,u,v] = loadEulerVelocities(fileName)

fid = fopen(fileName,'r');
val = fread(fid,'double');
fclose(fid);

nx = val(1);
ny = val(2);
val = val(3:end);

eulerX = zeros(nx,ny);
eulerY = zeros(nx,ny);
u = zeros(nx,ny);
v = zeros(nx,ny);

istart = 1;
% start of a pointer to where everything is stored in val


for k = 1:ny
  iend = istart + nx - 1;
  eulerX(1:nx,k) = val(istart:iend);
  istart = iend + 1;
end

for k = 1:ny
  iend = istart + nx - 1;
  eulerY(1:nx,k) = val(istart:iend);
  istart = iend + 1;
end

for k = 1:ny
  iend = istart + nx - 1;
  u(1:nx,k) = val(istart:iend);
  istart = iend + 1;
end

for k = 1:ny
  iend = istart + nx - 1;
  v(1:nx,k) = val(istart:iend);
  istart = iend + 1;
end

end % loadEulerVelocities

    
    