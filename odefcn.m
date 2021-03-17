function [dydt] = odefcn(t, params, N, rhatx, rhaty, rhatz, r, E, E0x, E0y, E0z, M, V, eps0, eps1, eps2, mu0, Lx, Ly, Lz, H0x, H0y, J, eta, R)
%This function represents the calculation of a differential equation that
%represents the change in angular velocity with respect ot time

D = zeros(3, N);
omega = zeros(3, N);

for i = 1:3*N
    if rem(i - 1, 3) == 0
        omega(1, floor(i/3) + 1) = params(i);
    elseif rem(i - 2, 3) == 0
        omega(2, floor(i/3) + 1) = params(i);
    else
        omega(3, floor(i/3)) = params(i);
    end
end

for i = (3*N + 1):6*N
    if rem(i - 1, 3) == 0
        D(1, floor((i - 3*N)/3) + 1) = params(i);
    elseif rem(i - 2, 3) == 0
        D(2, floor((i - 3*N)/3) + 1) = params(i);
    else
        D(3, floor((i - 3*N)/3)) = params(i);
    end
end

D = D./(ones(3,1)*sum(D.^2).^(1/2));

Dnx = ones(N,1)*D(1,:);
Dny = ones(N,1)*D(2,:);
Dnz = ones(N,1)*D(3,:);

Dmx = D(1,:)'*ones(1,N);
Dmy = D(2,:)'*ones(1,N);
Dmz = D(3,:)'*ones(1,N);

DnR = rhatx.*Dnx+rhaty.*Dny+rhatz.*Dnz;
DmR = rhatx.*Dmx+rhaty.*Dmy+rhatz.*Dmz;

DmDn = Dnx.*Dmx + Dny.*Dmy+ Dnz.*Dmz;

% ELECTRIC FIELD
% In this section, the local field directions are found by determining
% the local coordinates of each particle to write all local dipole
% directions, and transform back into global coordinates

% First use transformatins to determine local x,y,z coordinates
xold = [1,0,0];
yold = [0,1,0];
zold = [0,0,1];

znew = D(:,:,1);
xnew = zeros(N,3);
ynew = zeros(N,3);

rzq = zeros(N,4);

for iz = 1:N

    rzq(iz,:) = vrrotvec(zold,D(:,iz,1)); % find rotation axis and angle

    xnew(iz,:) = rodrigues_rot(xold,rzq(iz,1:3),rzq(iz,4)); % find local x-axis
    ynew(iz,:) = rodrigues_rot(yold,rzq(iz,1:3),rzq(iz,4)); % find local y-axis

end

% Next, dot product E0 with Dx, Dy, Dz for change in bases
E0nx = dot(xnew,E',2); % dotproduct of dx with E0  - Field in the local x-axis of the particle (coincides with a)
E0ny = dot(ynew,E',2); % Dot product of Dy with E0 -Field in the local y-axis of the particle (coincides with b)
E0nz = dot(znew,E)'; % Dot product of Dy with E0 - Field in the local z-axis of the particle (coincides with c)

plocx = V*(eps2 - eps1).*E0nx./(1 + ((eps2-eps1)/eps1)*Ly);
plocy = V*(eps2 - eps1).*E0ny./(1 + ((eps2-eps1)/eps1)*Ly);
plocz = V*(eps2 - eps1).*E0nz./(1 + ((eps2-eps1)/eps1)*Ly);

pglobx = plocx.*xnew(:,1) + plocy.*ynew(:,1) + plocz.*znew(1,:)';
pgloby = plocx.*xnew(:,2) + plocy.*ynew(:,2) + plocz.*znew(2,:)';
pglobz = plocx.*xnew(:,3) + plocy.*ynew(:,3) + plocz.*znew(3,:)';

pmR = rhatx.*pglobx + rhaty.*pgloby + rhatz.*pglobz;

Eddx = (1/4/pi/eps0)*sum((1./r.^3-eye(N)).*(3*pmR.*rhatx-pglobx)); % dipole-induced field effects in x-dir
Eddy = (1/4/pi/eps0)*sum((1./r.^3-eye(N)).*(3*pmR.*rhaty-pgloby)); % dipole-induced field effects in y-dir
Eddz = (1/4/pi/eps0)*sum((1./r.^3-eye(N)).*(3*pmR.*rhatz-pglobz)); % dipole-induced field effects in z-dir

% Electric torque calculations

% First calculate field at each particle location
Ex = E0x*ones(1,N) + Eddx; % total field in x
Ey = E0y*ones(1,N) + Eddy; % total field in y
Ez = E0z*ones(1,N) + Eddz; % total field in z

Eloc = [Ex; Ey; Ez];

% Transform field from global to local
Enx = dot(xnew,Eloc',2); % dotproduct of dx with E0  - Field in the local x-axis of the particle (coincides with a)
Eny = dot(ynew,Eloc',2); % Dot product of Dy with E0 -Field in the local y-axis of the particle (coincides with b)
Enz = dot(znew,Eloc)'; % Dot product of Dy with E0 - Field in the local z-axis of the particle (coincides with c)

% Determine local torques for each particle
Tex = V*((eps2 - eps1)^2).*(Lz - Ly).*Eny.*Enz./(eps1*(1 + ((eps2-eps1)/eps1)*Ly)*(1 + ((eps2-eps1)/eps1)*Lz));
Tey = V*((eps2 - eps1)^2).*(Lx - Lz).*Enx.*Enz./(eps1*(1 + ((eps2-eps1)/eps1)*Lx)*(1 + ((eps2-eps1)/eps1)*Lz));
Tez = V*((eps2 - eps1)^2).*(Ly - Lx).*Enx.*Eny./(eps1*(1 + ((eps2-eps1)/eps1)*Lx)*(1 + ((eps2-eps1)/eps1)*Ly));

% Global torques
Tnewx = Tex.*xnew(:,1) + Tey.*ynew(:,1) + Tez.*znew(1,:)';
Tnewy = Tex.*xnew(:,2) + Tey.*ynew(:,2) + Tez.*znew(2,:)';
Tnewz = Tex.*xnew(:,3) + Tey.*ynew(:,3) + Tez.*znew(3,:)';

% Set the magnetic field based on the current time
H0z = 10^5; %*sin(t*4*pi/10);

% Magnetic Field calculations
Hddx = (M*V/4/pi)*sum((1./r.^3-eye(N)).*(3*DmR.*rhatx-Dmx));
Hddy = (M*V/4/pi)*sum((1./r.^3-eye(N)).*(3*DmR.*rhaty-Dmy));
Hddz = (M*V/4/pi)*sum((1./r.^3-eye(N)).*(3*DmR.*rhatz-Dmz));

Hx = H0x*ones(1,N) + Hddx;
Hy = H0y*ones(1,N) + Hddy;
Hz = H0z*ones(1,N) + Hddz;

% Calculate magnetic torque on each particle
Tmx = mu0*M*V.*(Dny(1, :).*Hz - Dnz(1, :).*Hy);
Tmy = mu0*M*V.*(Dnz(1, :).*Hx - Dnx(1, :).*Hz);
Tmz = mu0*M*V.*(Dnx(1, :).*Hy - Dny(1, :).*Hx);

% Calculate the drag torque on each particle
Tdx = 8*pi*eta*R^3.*omega(1, :);
Tdy = 8*pi*eta*R^3.*omega(2, :);
Tdz = 8*pi*eta*R^3.*omega(3, :);

% Sum the torques in each direction 
Tx = Tmx - Tdx + Tnewx;
Ty = Tmy - Tdy + Tnewy;
Tz = Tmz - Tdz + Tnewz;

T = [Tx;Ty;Tz];

%dydt = zeros(6*N, 1);

omegaCrossD(1, :) = omega(2, :).*D(3, :) - omega(3, :).*D(2, :);
omegaCrossD(2, :) = omega(3, :).*D(1, :) - omega(1, :).*D(3, :);
omegaCrossD(3, :) = omega(1, :).*D(2, :) - omega(2, :).*D(1, :);

domega = omegaCrossD(1:3*N)';
dD = T(1:3*N)';

dydt = vertcat(domega, dD);

% for i = 1:6*N
%     if i <= 3*N
%         dydt(i) = omegaCrossD(i);
%     else
%         dydt(i) = T(i - 3*N)';
%     end
% end

end

