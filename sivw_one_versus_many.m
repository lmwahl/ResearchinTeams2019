% either vaccinate everyone against strain 10,
% or vaccinate uniformly against strains 8,9,10,11,12
% Compare results for these two cases, if epidemic starts at 10, 11, 12, ...
% Note that epidemics that start close to the boundaries of epidemic
% space are not treated properly, i.e. they will not get as big because
% we just drop any strains that mutate off the boundary

tic
clear totalI
nx = 20;                   % antigenic "classes" will run x = 1 to nx
deltat = .01;              % timestep
tend = 20;                 % final time
nsteps = tend/deltat;

plotflag = 1;              % whether to plot results each time

% parameter values, chosen pretty randomly
beta = 1;
sigma= 2;    % sd of mutation distribution
delta = 7;   % mean duration of infection, days 

ab = 2;   % half-sat for b function (infectibility reduced by vac)
hh = 2;    % hill coefficient for b
ad = 2;   % half-sat for d function (duration of infn reduced by vac)
jj = 2;    % hill coefficient for d
am = floor(nx/3);  % half-sat for m function (transmissibility changed)
kk = 2;    % hill coefficient for m
% using notation hh, jj and kk to avoid possible confusion with counters

dists =  0:1:nx;   % possible distances between two strains
% b function, increasing and saturating at beta
% annoyingly, matlab does not allow us to use 0 as an array index
% therefore bs(1) gives us b for antigenic distance 0 (likewise ds, mu, ms)
% and generally matlab bs(x) = value of function b evaluated at x-1
bs = beta*(dists).^hh./(dists.^hh + ab^hh);

% d function, currently starts at 1 and increases, saturating at delta
ds = (delta-1)*(dists).^jj./(dists.^jj + ad^jj);
ds = ds + 1;   %% bug fix because ds = 0 causes numerical issues %%

% mu as a function of antigenic distance, gaussian with sd sigma
% normalized to sum to 1
mu = exp(-(dists.*dists)./(sigma.^2));
mu = mu./(mu(1) +2*sum(mu(2:end)));

% m function, mu of transmitted strains modified by vaccination
for x=1:nx
  for y = 1:nx
    for z = 1:nx
       xz = abs(x-z);
       ms(x,y,z) = mu(abs(x-y)+1)*(xz^kk/(xz^kk + am^kk));
    end
  end
end

% initializing arrays to the correct size
Sinit = 0;
Iinit = zeros(1,nx); 
Winit = zeros(1,nx,nx);
Vinit = zeros(1,nx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial conditions (edit here)

for iinit = 10:20
for caseflag = 1:2
 S = Sinit; I = Iinit; W = Winit; V = Vinit;
 I(iinit) = .01;
 if (caseflag==1)
  V(10) = .99;
 else
  V(8:12) = .99/5;  
 end

for istep = 2:nsteps    % loop for the numerical integration

  %%%%%%%%%%%  dS/dt %%%%%%%%%%%%%%
  Isum = sum(I(istep-1,:));
  Wsum = 0;
  for iy = 1:nx
  for iz = 1:nx
    msum = 0;
    for ix = 1:nx, msum = msum + ms(ix,iy,iz); end
    Wsum = Wsum + W(istep-1,iy,iz)*msum;
  end
  end
  S(istep) = S(istep-1) - beta*S(istep-1)*(Isum + Wsum)*deltat;

  %%%%   dI/dt %%%%%%%%%%%%%%
  for ix = 1:nx
    Imusum = 0;
    for iy =1:nx
      Imusum = Imusum + I(istep-1,iy).*mu(abs(iy-ix)+1);
    end
    mWsum = 0;
    for iz = 1:nx
    for iy = 1:nx
      mWsum = mWsum + ms(ix,iy,iz)*W(istep-1,iy,iz);
    end
    end
    I(istep,ix) = I(istep-1,ix) + beta*S(istep-1)*(Imusum + mWsum)*deltat;
    I(istep,ix) = I(istep,ix) - I(istep-1,ix)/delta*deltat;
  end

  %%%%%%%%%%%% dV/dt %%%%%%%%%%%%%%%%%
  for iz = 1:nx
    bImusum = 0;
    bmWsum = 0;
    for iy = 1:nx
     for ix = 1:nx
        bImusum = bImusum + I(istep-1,iy)*bs(abs(iz-ix)+1)*mu(abs(ix-iy)+1);
      for iu = 1:nx
        bmWsum = bmWsum + W(istep-1,iy,iu)*bs(abs(iz-ix)+1)*ms(ix,iy,iu);
      end
     end
    end
    V(istep,iz) = V(istep-1,iz) - V(istep-1,iz)*(bImusum + bmWsum)*deltat;
  end

  %%%%%%%%%%  dW/dt %%%%%%%%%%%%%%%%%%%%%%%5
  for ix = 1:nx
  for iz = 1:nx
    Imusum = 0;
    for iy =1:nx
      Imusum = Imusum + I(istep-1,iy)*mu(abs(iy-ix)+1);
    end
    mWsum = 0;
    for iu = 1:nx
    for iy = 1:nx
      mWsum = mWsum + ms(ix,iy,iu)*W(istep-1,iy,iu);
    end
    end
    W(istep,ix,iz) = W(istep-1,ix,iz) + bs(abs(iz-ix)+1)*V(istep-1,iz)*(Imusum + mWsum)*deltat;
    W(istep,ix,iz) = W(istep,ix,iz) - W(istep-1,ix,iz)/ds(abs(ix-iz)+1)*deltat;
  end
  end
  end  %%% loop on istep

% plot some graphs (turn off for many reps)
if (plotflag)
 Winfectedt = sum(W,3);
 Wvaccinet  = sum(W,2);
 Wvaccinet = reshape(Wvaccinet,nsteps,nx);
 figure(1);  
 mesh(I)
 figure(2);
 mesh(Winfectedt);
 title('W infected by');
 figure(3);
 mesh(Wvaccinet);
 title('W vaccinated against');
end

% compute total number of infection-days
totalI(iinit,caseflag) = sum(sum(sum(W))) + sum(sum(I));

end % loop on caseflag
end % loop on iinit

plot(totalI,'-o');
legend('vaccine concentrated at 10','vaccine distributed 8-12','location','NW');
xlabel('incoming strain');
ylabel('total infection days');
set(gca,'XLim',[10 18])
toc
