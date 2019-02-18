% do multiple simulations per initial location

absorbing_radius = 2; % radius of the inner(absorbing) circle
bounding_radius = 4; % raduis of the outer(boundary) circle


r = @(theta) absorbing_radius; %our function for our absorbing surface (for plotting purposes)
r2 = @(theta) bounding_radius; %our function for our boundary circle (for plotting purposes)

% this block of code is for plotting our circles
theta = linspace(0,2*pi,1000);
geom = r(theta).*exp(1i*theta);
geom2 = r2(theta).*exp(1i*theta);



nx = 50; %number of points in the x direction
ny = 50; %number of points in the y direction

%!!!DO NOT MAKE ANY OF THESE NUMBERS LESS THAN THE BOUNDARY RADIUS!!!
xmin = -5;
xmax = 5;
ymin = -5;
ymax = 5;

[x,y] = meshgrid(linspace(xmin,xmax,nx),linspace(ymin,ymax,ny));

init_pts = x(:) + 1i*y(:);% initial location of the random walks
npts = numel(init_pts);% number of initial points

nruns = 100;% number of runs per initial location
mfpt = zeros(npts,1);% initialize our mean first passage time vector

D = 1; % diffusion coefficient

dt = 0.1;% our step size for the random walker
for k = 1:numel(init_pts)%iterate over each of the nx by ny inital pts
      if abs(init_pts(k)) > bounding_radius %if the initial location of the random walker...
                                            %is not in the boundary circle then...  
          % ...Do nothing
      else    
      disp(k)
      z = init_pts(k)*ones(nruns,1); 
      zold = z;
      znew = zold;
      dist = sqrt(real(z).^2+imag(z).^2); %since our absorbing circle is the circle with center at the origin
      fpt = zeros(nruns,1);
      while any(dist > absorbing_radius) 
        s = find(dist > absorbing_radius);
        fpt(s) = fpt(s) + dt;
        omega = 2*pi*rand(nruns,1);
        znew(s) = zold(s) + D*dt*exp(1i*omega(s));
    
        if any(abs(znew) > bounding_radius) % to account for our closed boundary
          s = find(abs(znew) > bounding_radius);
          overshoot = abs(znew(s)) - bounding_radius;
          znew(s) = (bounding_radius - overshoot).*znew(s)./abs(znew(s));
        end
        
        dist = sqrt(real(znew).^2+imag(znew).^2);
        zold = znew;
      end
      
      mfpt(k) = mean(fpt); % stores the mean first passage of the random walker
      end
end

mfpt = reshape(mfpt,ny,nx);

clf;
surf(x,y,mfpt);
colorbar
view(2); shading interp;
axis equal; hold on
plot3(real(geom),imag(geom),1000*ones(size(geom)),'k','linewidth',2)
plot3(real(geom2),imag(geom2),1000*ones(size(geom2)),'k','linewidth',2)
grid on
 
axis([xmin xmax ymin ymax])
hold off

