% do multiple simulations per initial location

% poisonville
vec1 = @(z) 1 - imag(z).^2;

% extensional
vec2 = @(z) -real(z) + 1i*imag(z);

% counterclockwise rotation
vec3 = @(z) 1i*z;

% shear
vec4 = @(z) imag(z);

absorbing_radius = 1/3;
N = 5;

%!!!DO NOT MAKE ANY OF THESE NUMBERS LESS THAN THE BOUNDARY RADIUS!!!
xmin = -1;
xmax = 1;
ymin = -1;
ymax = 1;

bounding_radius = 1; % raduis of the outer(boundary) circle
[xc, yc] = center_generator_rect(absorbing_radius,N,xmin,xmax,ymin,ymax);
centers = xc + 1i*yc;

% this block of code is for plotting our circles
theta = linspace(0,2*pi,1000);

geoms = zeros(1000, N);
for i=1:N
    disp( [xc(i) yc(i)] )
    geoms(:, i) = absorbing_radius*exp(1i*theta) + centers(i);    
end

nx = 50; %number of points in the x direction
ny = 50; %number of points in the y direction


[x,y] = meshgrid(linspace(xmin,xmax,nx),linspace(ymin,ymax,ny));

init_pts = x(:) + 1i*y(:);% initial location of the random walks
npts = numel(init_pts);% number of initial points

nruns = 50;% number of runs per initial location
mfpt = zeros(npts,1);% initialize our mean first passage time vector
 
D = 1; % diffusion coefficient
 
dt = 0.01;% our step size for the random walker
for k = 1:numel(init_pts)%iterate over each of the nx by ny inital pts  
      disp(init_pts(k))
      z = init_pts(k)*ones(nruns,1); 
      zold = z;
      znew = zold;
      
      dist = abs(z - centers(1));
      for i=2:N
          dist = min( dist, abs(z - centers(i)));
      end
      
      fpt = zeros(nruns,1);
      % min_dist.m
      while any(dist > absorbing_radius) 
        s = find(dist > absorbing_radius);
        
        fpt(s) = fpt(s) + dt;
        omega = 2*pi*rand(nruns,1);
        znew(s) = zold(s) + D*dt*(exp(1i*omega(s)) + vec3(zold(s)));
    
        if any((imag(znew) > ymax) + (imag(znew) < ymin) + ...
               (real(znew) > xmax) + (real(znew) < ymin) ) % to account for our closed boundary
            s = find(abs(znew) > bounding_radius);    
            % Simulates the bounding radius acting as an absorbing radius
            znew(s) = centers(1);
            % overshoot = abs(znew(s)) - bounding_radius;
            % znew(s) = (bounding_radius - overshoot).*znew(s)./abs(znew(s));
        end
        % end min_dist
        dist = abs(znew - centers(1));
        for i=2:N
          dist = min( dist, abs(znew - centers(i)));
        end

        zold = znew;
      end
      
      mfpt(k) = mean(fpt); % stores the mean first passage of the random walker
      
end

ss = find(mfpt < 0.01*dt);
mfpt(ss) = 1/0;
mfpt = reshape(mfpt,ny,nx);

clf;
surf(x,y,mfpt);
colorbar
view(2); shading interp;
axis equal; hold on
for i=1:N
    fill(real(geoms(:,i)),imag(geoms(:,i)),'k');
end

grid off
axis([xmin xmax ymin ymax])
set(gca,'visible', 'off'); 
hold off
