nx = 20;
deltat = .01;
tend = 20;
nsteps = tend/deltat;

nreps = 25;
plotflag = 0;

beta = 1;
sigmamu = 2;
delta = 7;
af = floor(nx/3);
ag = 3;
ah = floor(nx/4);
xs = 0:1:nx;
qs = exp(xs)./(exp(xs) + exp(ah));
%qs = ones(size(xs));  
ds = delta*exp(xs)./(exp(xs) + exp(af));
ds = ds + 1;   %%%%%%%%%%fix this later %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
mu = exp(-(xs.*xs)./(sigmamu.^2));
mu = mu./(sum(mu));
for x=1:nx
  for y = 1:nx
    for z = 1:nx
      alpha(x,y,z) = mu(abs(x-y)+1)*(exp(abs(x-z))/(exp(abs(x-z)) + exp(ag)));
    end
  end
end

prevcentre = 10;
newcentre= 12;

for irep = 1:nreps

S = 0;
%S = .99;
I = zeros(1,nx);
for iinits = 1:5
  istrain(irep,iinits) = round(10 + sigmamu*randn);
  I(istrain(irep,iinits)) = I(istrain(irep,iinits)) + .01/5;
end
V = zeros(1,nx);
V(8:12) = .99/5;
%V(10) = .99;
W = zeros(1,nx,nx);

for istep = 2:nsteps
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

totalI(irep) = sum(sum(sum(W))) + sum(sum(I));

end  % loop on irep

istrain
totalI
mean(totalI)
std(totalI)
