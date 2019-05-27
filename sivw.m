tic
%%% this needs a re-write, currently in the "old" notation

nx = 20;                   % antigenic "classes" will run x = 1 to nx
deltat = .01;              % timestep
tend = 20;                 % final time
nsteps = tend/deltat;

nreps = 25;                % number of stochastic runs
ninits = 5;                % number of independent introductions that will
                           % start the epidemic this season
plotflag = 0;              % whether to plot results each time

% parameter values, chosen pretty randomly
beta = 1;
sigmamu = 2;
delta = 7;
af = floor(nx/3);
ag = 3;
ah = floor(nx/4);

xs = 0:1:nx;
% q function, increasing and saturating at 1
% annoyingly, matlab does not allow us to use 0 as an array index
% therefore qs(1) = q for antigenic distance 0 (likewise ds, mu, alpha)
% and generally matlab qs(x) = value of function q evaluated at x-1
qs = exp(xs)./(exp(xs) + exp(ah));

% d function, currently starts at 1 and increases, saturating at delta
ds = delta*exp(xs)./(exp(xs) + exp(af));
ds = ds + 1;   %% bug fix because ds = 0 causes numerical issues %%

% mu as a function of antigenic distance, gaussian with sd sigmamu
% normalized to sum to 1
mu = exp(-(xs.*xs)./(sigmamu.^2));
mu = mu./(sum(mu));

% alpha function, see the manuscript for explanation!
for x=1:nx
  for y = 1:nx
    for z = 1:nx
      alpha(x,y,z) = mu(abs(x-y)+1)*(exp(abs(x-z))/(exp(abs(x-z)) + exp(ag)));
    end
  end
end

% when simulating a new season, with the incoming transmission cloud
% centred differently than the previous season
prevcentre = 8;   % most common strain last season

% set up "saved" initial conditions to reuse in each rep
Sinit = 0;    % no one is initially susceptible
Iinit = zeros(1,nx);  % no one is initially infected
Winit = zeros(1,nx,nx);
Vinit = zeros(1,nx);

%  99% of population is vaccinated.  Use 5 antigenic positions centred around
% prevcentre
%
%Vinit(prevcentre-2:prevcentre+2) = .99/5;
%fprintf(1,'Distributed vaccination, centred at strain from last year\n');

%  OR commented out below, just vaccinate everyone at prevcentre
%
Vinit(prevcentre) = .99;
fprintf(1,'Vaccinating everyone with strain from last year\n');


for irep = 1:nreps    %loop over different stochastic trials

% choose where the epidemic this year is centred
% currently doubling the s.d. of the mu distribution to get more
% variance to try to see an effect
newcentre= round(prevcentre + 2*sigmamu*randn);
newcentre = max(newcentre,1); newcentre = min(newcentre,nx);
newcs(irep) = newcentre;  % save these to check later

S = Sinit; I = Iinit; W = Winit; V = Vinit;

%ninits new strains are going to simultaneously
% start the season off.  Chosen centred around newcentre 
for iinits = 1:ninits   
  % what new strains will rain down from the cloud?  store them in istrain
  tempstrain = round(newcentre + sigmamu*randn);
  tempstrain = max(tempstrain,1); tempstrain = min(tempstrain,nx);
  istrain(irep,iinits) = tempstrain;
  % in total 1% of the population will be infected, 1/ninits of that for each strain
  I(istrain(irep,iinits)) = I(istrain(irep,iinits)) + .01/ninits;
end


for istep = 2:nsteps    % loop for the numerical integration
  %%%%%%%%%%%  dS/dt %%%%%%%%%%%%%%
  Isum = sum(I(istep-1,:));
  Wsum = 0;
  for iy = 1:nx
  for iz = 1:nx
    alsum = 0;
    for ix = 1:nx, alsum = alsum + alpha(ix,iy,iz); end
    Wsum = Wsum + W(istep-1,iy,iz)*alsum;
  end
  end
  S(istep) = S(istep-1) - beta*S(istep-1)*(Isum + Wsum)*deltat;

  %%%%   dI/dt %%%%%%%%%%%%%%
  for ix = 1:nx
    Imusum = 0;
    for iy =1:nx
      Imusum = Imusum + I(istep-1,iy).*mu(abs(iy-ix)+1);
    end
    alphaWsum = 0;
    for iz = 1:nx
    for iy = 1:nx
      alphaWsum = alphaWsum * alpha(ix,iy,iz)*W(istep-1,iy,iz);
    end
    end
    I(istep,ix) = I(istep-1,ix) + beta*S(istep-1)*(Imusum + alphaWsum)*deltat;
    I(istep,ix) = I(istep,ix) - I(istep-1,ix)/delta*deltat;
  end

  %%%%%%%%%%%% dV/dt %%%%%%%%%%%%%%%%%
  for iz = 1:nx
    qImusum = 0;
    qWalphasum = 0;
    for iy = 1:nx
     for ix = 1:nx
        qImusum = qImusum + I(istep-1,iy)*qs(abs(iz-ix)+1)*mu(abs(ix-iy)+1);
      for iu = 1:nx
        qWalphasum = qWalphasum + W(istep-1,iy,iu)*qs(abs(iz-ix)+1)*alpha(ix,iy,iu);
      end
     end
    end
    V(istep,iz) = V(istep-1,iz) - beta*V(istep-1,iz)*(qImusum + qWalphasum)*deltat;
  end

  %%%%%%%%%%  dW/dt %%%%%%%%%%%%%%%%%%%%%%%5
  for ix = 1:nx
  for iz = 1:nx
    Imuqsum = 0;
    for iy =1:nx
      Imuqsum = Imuqsum + I(istep-1,iy)*mu(abs(iy-ix)+1);
    end
    alphaWqsum = 0;
    for iu = 1:nx
    for iy = 1:nx
      alphaWqsum = alphaWqsum * alpha(ix,iy,iu)*W(istep-1,iy,iu)*qs(abs(iz-ix)+1);
    end
    end
    W(istep,ix,iz) = W(istep-1,ix,iz) + beta*qs(abs(iz-ix)+1)*V(istep-1,iz)*(Imuqsum + alphaWqsum)*deltat;
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
totalI(irep) = sum(sum(sum(W))) + sum(sum(I));

end  % loop on irep

% print stuff out to look at, or not
newcs
%istrain
%totalI
fprintf('Mean number of infection days %f +/- %f\n',mean(totalI),std(totalI))
toc
