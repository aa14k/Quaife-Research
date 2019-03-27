% do multiple simulations per initial location

absorbing_radius = 1;
N = 3;

%!!!DO NOT MAKE ANY OF THESE NUMBERS LESS THAN THE BOUNDARY RADIUS!!!
xmin = -3;
xmax = 3;
ymin = -3;
ymax = 3;

bounding_radius = 3; % raduis of the outer(boundary) circle
[xc, yc] = center_generator(absorbing_radius,N,xmin,xmax,ymin,ymax,bounding_radius);
centers = xc + 1i*yc;

%r = @(theta) absorbing_radius; %our function for our absorbing surface (for plotting purposes)
r2 = @(theta) bounding_radius; %our function for our boundary circle (for plotting purposes)

% this block of code is for plotting our circles
theta = linspace(0,2*pi,1000);

geoms = zeros(1000, N);
for i=1:N
    disp( [xc(i) yc(i)] )
    geoms(:, i) = exp(1i*theta) + centers(i);    
end
geom2 = r2(theta).*exp(1i*theta);



nx = 400; %number of points in the x direction
ny = 400; %number of points in the y direction


[x,y] = meshgrid(linspace(xmin,xmax,nx),linspace(ymin,ymax,ny));

init_pts = x(:) + 1i*y(:);% initial location of the random walks
npts = numel(init_pts);% number of initial points

nruns = 50;% number of runs per initial location
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
      
      dist = abs(z - centers(1));
      for i=2:N
          dist = min( dist, abs(z - centers(i)));
      end
          %dist = min(abs(z-centers(1, 1)), abs(z-centers(2, 1)), abs(z-centers(3, 1)));%sqrt(real(z).^2+imag(z).^2); %since our absorbing circle is the circle with center at the origin
      
      %[inds, dist] = dsearchn(centers, z);
      %dist = abs(dist);
      
      fpt = zeros(nruns,1);
      % min_dist.m
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
        % end min_dist
        dist = abs(znew - centers(1));
        for i=2:N
          dist = min( dist, abs(znew - centers(i)));
        end
        %dist = min(abs(znew-centers(1, 1)), abs(znew-centers(2, 1)), abs(znew-centers(3, 1)));
        %[inds, dist] = dsearchn(centers, znew);
        %dist = abs(dist);
        zold = znew;
      end
      
      mfpt(k) = mean(fpt); % stores the mean first passage of the random walker
      end
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
plot3(real(geom2),imag(geom2),1000*ones(size(geom2)),'k','linewidth',2, 'Color', 'k')

grid off
axis([xmin xmax ymin ymax])
set(gca,'visible', 'off'); 
hold off
