% do multiple simulations per initial location

absorbing_radius = 2;
bounding_radius = 2*absorbing_radius;

% all these are to plot different circles for perspective
r = @(theta) absorbing_radius;
r2 = @(theta) 2*absorbing_radius;
r3 = @(theta) 3*absorbing_radius;
r4 = @(theta) 4*absorbing_radius;
r5 = @(theta) 5*absorbing_radius;

theta = linspace(0,2*pi,1000);
geom = r(theta).*exp(1i*theta);
geom2 = r2(theta).*exp(1i*theta);
geom3 = r3(theta).*exp(1i*theta);
geom4 = r4(theta).*exp(1i*theta);
geom5 = r5(theta).*exp(1i*theta);

%number of initial points on the xy plane
nx = 100;
ny = 100;

%the boundary of the space, do not make any of these values less than the bounding_radius.
xmin = -10;
xmax = 10;
ymin = -10;
ymax = 10;
[x,y] = meshgrid(linspace(xmin,xmax,nx),linspace(ymin,ymax,ny));

init_pts = x(:) + 1i*y(:);
% initial location of the random walks
npts = numel(init_pts);
% number of initial points

nruns = 100;
% number of runs per initial location
mfpt = zeros(npts,1);

D = 1; % diffusion coefficient


dt = 0.1;
for k = 1:numel(init_pts)
  disp(k)
  z = init_pts(k)*ones(nruns,1);
  zold = z;
  znew = zold;

  dist = sqrt(real(z).^2+imag(z).^2); %since our geometry is a circle at the origin
  fpt = zeros(nruns,1);
  while any(dist > absorbing_radius)
    s = find(dist > absorbing_radius);
%    disp(numel(s))
    fpt(s) = fpt(s) + dt;
    omega = 2*pi*rand(nruns,1);
    znew(s) = zold(s) + D*dt*exp(1i*omega(s));

    if any(abs(znew) > bounding_radius)
      s = find(abs(znew) > bounding_radius);
      overshoot = abs(znew(s)) - bounding_radius;
      znew(s) = (bounding_radius - overshoot).*znew(s)./abs(znew(s));
    end

    dist = sqrt(real(znew).^2+imag(znew).^2); %since our geometry is a circle at the origin

    zold = znew;
  end
  mfpt(k) = mean(fpt);
end

mfpt = reshape(mfpt,ny,nx);

clf;
surf(x,y,mfpt);
view(2); shading interp;
axis equal; hold on
plot3(real(geom),imag(geom),1000*ones(size(geom)),'k','linewidth',2)
plot3(real(geom2),imag(geom2),1000*ones(size(geom2)),'k','linewidth',2)
plot3(real(geom3),imag(geom3),1000*ones(size(geom3)),'k','linewidth',2)
plot3(real(geom4),imag(geom4),1000*ones(size(geom4)),'k','linewidth',2)
plot3(real(geom5),imag(geom5),1000*ones(size(geom5)),'k','linewidth',2)

axis([xmin xmax ymin ymax])


