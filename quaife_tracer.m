%% Loading Data
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

% Create interpolating objects that can be queried at given points
Fx = griddedInterpolant( eulerX', eulerY', u', 'linear' );
Fy = griddedInterpolant( eulerX', eulerY', v', 'linear' );

%% Bounding Shape
xmin = 2;
xmax = 3;
ymin = 2;
ymax = 3;

% Plot only a subsection of the circles
sL = real(centers) > xmin - 0.3;
sR = real(centers) < xmax + 0.3;
sU = imag(centers) > ymin - 0.3;
sD = imag(centers) < ymax + 0.3;

centers = centers( sL & sR & sU & sD );
radii = radii( sL & sR & sU & sD );

N = numel(centers);

%% Parameters   
p1 = 1; % Diffusion
p2 = 1; % Advection
dt = 0.1;% our step size for the random walker

%% Tracers
nruns = 1000;% number of runs per initial location

% Initialize tracer runs
%z = [1.2+1.8i ; 1.4+1.2i; 1.3+1.5i; 1.6+1.2i];
%z = (1.2 + 1.8i)*ones(nruns,1); 
z = xmin + 1i*linspace(ymin, ymax, nruns)';
z = repmat( z, 3, 1 );

%z = 0.6 + 1.4i;

%th = 0:pi/500:2*pi;
%xunit = 0.2 * cos(th) + 1.6;
%yunit = 0.2 * sin(th) + 1.4;
%z = xunit' + 1i*yunit';

nruns = numel( z );
fpt = zeros(nruns,1); % vector of first passage times
count_absorbed = 0;   % Number of tracers that hit other edge

z_paths = z(:);
zold = z;
znew = zold;

%% Check which are inside circles
inside = zeros(size(z));
for i=1:nruns
    inside(i) = any( abs( znew(i) - centers ) < radii );
end

while any(~inside) 
    s = find(~inside);
    disp( sum( inside ) );
    
    fpt(s) = fpt(s) + dt;
    %omega = 2*pi*rand(nruns,1);
    %alpha = rand(nruns, 1);
    gamma = normrnd(0, sqrt(dt), nruns, 2);
    
    xinterp = Fx( real(zold(s)), imag(zold(s)) );
    yinterp = Fy( real(zold(s)), imag(zold(s)) );

    znew(s) = zold(s) + p1*dt*(gamma(:,1) + 1i*gamma(:,2)) + ...
                        p2*dt*(xinterp + 1i*yinterp);
    
    %znew(s) = zold(s) + p1*dt*exp(1i*omega(s)) .* sqrt(-dt*2*log(alpha(s))) + ...
    %                    p2*dt*(xinterp + 1i*yinterp);

    z_paths = [z_paths znew];

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
        count_absorbed = count_absorbed + 1;
    end


    zold = znew;
end


disp( 'mfpt:' )
disp(mean(fpt)); % stores the mean first passage of the random walker
disp( 'absorbed by edge' )
disp( count_absorbed )

% plot circles
figure(2);
clf;
axis equal; hold on

for i=1:nruns
    if( real(z_paths(i,end)) >= xmax || ...
        real(z_paths(i,end)) <= xmin || ...
        imag(z_paths(i,end)) >= ymax || ...
        imag(z_paths(i,end)) <= xmin)
        p1 = plot( real(z_paths(i,:)), imag(z_paths(i,:)), 'b-' );        
    else
        p1 = plot( real(z_paths(i,:)), imag(z_paths(i,:)), 'r-' );
    end
    p1.Color(4) = 0.2;
end

scatter( real(z), imag(z), 'r.' );

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

    
    