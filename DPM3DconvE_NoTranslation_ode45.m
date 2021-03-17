function [data] = DPM3DconvE_NoTranslation_ode45(defaultN, volFrac, varargin)
%% Input Parsing

% This function was created to accept different types of input:
tic
p = inputParser;

% Any input that is not defined will take its corresponding following
% value:
defaultN    = defaultN;%5 ; % numberof particles
defaultvol  = volFrac;
defaultHtoM = 5; % ratio of magnetic to electric fields
defaultR    = 0.5*10^-6;
defaultrhok = 5.28*(1/1000)*(100/1)^3;% density of particles
defaultM    = 72*defaultrhok; % magnetization
defaulteta  = 3500*(1/1000);%dieletric permitivitty of XX
defaulttmax = 6000000*3;
defaultmaxit= 10000;

% This function can take as input, its own output. If this option is used,
% the simulation will continue wherever it left off the previous run.
defaultdata = 0;

addOptional(p,'N',defaultN,@isnumeric);
addOptional(p,'vol',defaultvol,@isnumeric);
addOptional(p,'R',defaultR,@isnumeric);
addOptional(p,'M',defaultM,@isnumeric);
addOptional(p,'HtoM',defaultHtoM,@isnumeric);
addOptional(p,'eta',defaulteta,@isnumeric);
addOptional(p,'tmax',defaulttmax,@isnumeric);
addOptional(p,'maxit',defaultmaxit,@isnumeric);
addOptional(p,'data',defaultdata);

parse(p,varargin{:});

%% Initialization

datain = p.Results.data;
tmax = p.Results.tmax;
maxit= p.Results.maxit;

if isnumeric(datain)==0
    
    Param = datain.Parameters;
    N     = Param(1);
    vol   = Param(2);
    HtoM  = Param(3);
    R     = Param(4);
    dim   = Param(5);
    M     = Param(6);
    eta   = Param(7);
    
    t = datain.t1(end);
    
    X = datain.X(:,:,end);
    D = datain.D(:,:,end);
    
    Xdata       = datain.X;
    Ddata       = datain.D;
    tdata       = datain.t1;
    g           = datain.g;
    f           = datain.f;
    previous_i  = length(datain.t1);
  
else
    
    N    = p.Results.N;
    vol  = p.Results.vol;
    HtoM = p.Results.HtoM;
    R    = p.Results.R;
    dim  = (N*(4/3)*pi*R^3/vol)^(1/3)/2;
    M    = p.Results.M;
    eta  = p.Results.eta;

    t      = 0;
    
    fid = fopen('initData.txt','wt');
    fprintf(fid,'%g\t %g\t %g\t' ,N,R,dim);
    fprintf(fid,'\n');
    fclose(fid);
    
    x=randomhardparticles3D(N,R,dim);
    %X(:,:)    = permute(x(:,1,:),[1,3,2]);
    D(:,:)    = permute(x(:,2,:),[1,3,2]);
%     d1 = ones(2, N)*sqrt(0.5);
%     d2 = zeros(1, N);
%     D = vertcat(d1, d2);
    
    X = ones(3, N); 
    interval = ceil((N)^(1/3));
    step = 2*dim/(interval);
    index = 1;
    for i = 1:interval
        for j = 1:interval
            for k = 1: interval
                if index <= N
                    %X(:, index) = [-dim + step*k, -dim + step*j, -dim + step*i];
                    X(:, index) = [-dim + step/2 + step*(k-1), -dim + step/2 + step*(j-1), -dim + step/2 + step*(i-1)];
                    index = index + 1;
                end
            end
        end
    end
    
    saveFile='test';
    resPos = strcat(saveFile,'_initPos.txt');
    fid = fopen(resPos,'wt');
    for ii = 1:size(X,1)
        fprintf(fid,'%g\t',X(ii,:));
        fprintf(fid,'\n');
    end
    fclose(fid);
    
    resDir = strcat(saveFile,'_initDir.txt');
    fid = fopen(resDir,'wt');
    for ii = 1:size(D,1)
        fprintf(fid,'%g\t',D(ii,:));
        fprintf(fid,'\n');
    end
    fclose(fid);

    
    Xdata(:,:,1)  = X;
    Ddata(:,:,1)  = D;
    g(1)          = zeros(1);
    f(1)          = zeros(1);
    tdata(1)      = 0;
    previous_i    = 1;
end

% Plotting conditions are defined
plotmult = 1000; % this determines the multiples of it at which the program will plot the state
times = zeros(maxit/plotmult,1);

%%%%%%%% APPLY ELELCTRIC AND MAGNETIC EXTERNAL FIELDS %%%%%%%%%%%
E0x = 0; E0y = 0; E0z = 0; %100*10^6; %applied ELECTRIC field [V/m]
H0x = 0; H0y = 0; H0z = 0; %HtoM*M; %applied MAGNETIC field [A/m]

% Declare initial condition for angular velocity
omega = zeros(3, N);

mu0   = 4*pi*10^(-7); %permeability of free space [N/A^2]
B     = mu0*H0z; %induced field

%E = [E0x, E0y, E0z]; % initial field directions
% Initial field setup
    Ex = E0x*ones(1,N); % total field in x
    Ey = E0y*ones(1,N); % total field in y
    Ez = E0z*ones(1,N); % total field in z
    E = [Ex; Ey; Ez];

%Particle and magnetic field properties:
Vequi    = 4/3*pi*R^3; %particle volume [m^3] 
a = ((Vequi*(3/4)*3.5)/(pi))^(1/3);
b = a;
c = a/3.5;

V    = 4/3*pi*a*b*c; % [m^3] ellipsoidal particle volume

Lx = (a^2*c)/(2*(a^2 - c^2))*(((pi/2)/sqrt(a^2-c^2))-(c/(a^2))); % depolarizing factor in x-dir
Ly = Lx; % depolarizing factor in y-dir
Lz = 1 - 2*Lx; % depolarizing factor in z-dir

eps0 = 8.85*10^-12; % [F/m] electric permittivity of vacuum
eps1 = 2*eps0; % electric permittivity of medium
eps2 = 10*eps0; % electric permittivity of particles

% Magnetic field
%Ds    = 6*pi*eta*R; %Drag coefficient
eta   = p.Results.eta; %viscosity of the fluid [Pa*s]

Q    = mu0*(norm(M)*V)^2/(4*pi); %Magnetic force constant
QQQ  = mu0*M*V/(8*pi*eta*R^3); % Magnetic Torque constant
QQ   = 3*Q*norm(mu0*HtoM*M)^2/(8*R^4); %Repelling force constant
QQQe = 1/(8*pi*eta*R^3); % Electric torque
beta = 40; %Repelling parameter
J(1) = 2/5*R^3; %Moment of inertia of particles

% Useful operators
corrector  = ones(N,N)-2*tril(ones(N,N),-1);
corrector2 = diag(inf(1,N));

dim2 = dim;
%diffx=@(XI,XJ) diag(([XI(1)-XJ(:,1),XI(1)-XJ(:,1)+2*dim,XI(1)-XJ(:,1)-2*dim])*(abs([XI(1)-XJ(:,1),XI(1)-XJ(:,1)+2*dim,XI(1)-XJ(:,1)-2*dim])<dim2)');
%diffy=@(XI,XJ) diag(([XI(2)-XJ(:,2),XI(2)-XJ(:,2)+2*dim,XI(2)-XJ(:,2)-2*dim])*(abs([XI(2)-XJ(:,2),XI(2)-XJ(:,2)+2*dim,XI(2)-XJ(:,2)-2*dim])<dim2)');
%diffz=@(XI,XJ) diag(([XI(3)-XJ(:,3),XI(3)-XJ(:,3)+2*dim,XI(3)-XJ(:,3)-2*dim])*(abs([XI(3)-XJ(:,3),XI(3)-XJ(:,3)+2*dim,XI(3)-XJ(:,3)-2*dim])<dim2)');
diffx=@(XI,XJ) diag(([XI(1)-XJ(:,1),XI(1)-XJ(:,1)+2*dim,XI(1)-XJ(:,1)-2*dim])*(([abs(XI(1)-XJ(:,1))<=dim2,abs(XI(1)-XJ(:,1)+2*dim)<dim2,abs(XI(1)-XJ(:,1)-2*dim)<dim2])'));
diffy=@(XI,XJ) diag(([XI(2)-XJ(:,2),XI(2)-XJ(:,2)+2*dim,XI(2)-XJ(:,2)-2*dim])*(([abs(XI(2)-XJ(:,2))<=dim2,abs(XI(2)-XJ(:,2)+2*dim)<dim2,abs(XI(2)-XJ(:,2)-2*dim)<dim2])'));
diffz=@(XI,XJ) diag(([XI(3)-XJ(:,3),XI(3)-XJ(:,3)+2*dim,XI(3)-XJ(:,3)-2*dim])*(([abs(XI(3)-XJ(:,3))<=dim2,abs(XI(3)-XJ(:,3)+2*dim)<dim2,abs(XI(3)-XJ(:,3)-2*dim)<dim2])'));

xunit=[ones(1,N);zeros(1,N);zeros(1,N)];
yunit=[zeros(1,N);ones(1,N);zeros(1,N)];
zunit=[zeros(1,N);zeros(1,N);ones(1,N)];

i = previous_i+1;
it= 1;
tf= 0;
ti= cputime;
h = 0; % ADDED - REMOVE LATER
maxDeflection = 0;

%% Iterative Solver
%while it<maxit && tf<tmax*60 && t<0.5
    
%while t<i*0.0001 && t<0.5
    
%while it<maxit && tf<tmax*60 && maxDeflection < 10
    
%while t<i*0.0001 && maxDeflection < 10

    %Boundary condition implementation

    Ixplus=find(dot(X,xunit)>dim);
    Ixminus=find(dot(X,xunit)<-dim);
    Iyplus=find(dot(X,yunit)>dim);
    Iyminus=find(dot(X,yunit)<-dim);
    Izplus=find(dot(X,zunit)>dim);
    Izminus=find(dot(X,zunit)<-dim);
        
    X(1,Ixplus)=X(1,Ixplus)-2*dim;
    X(1,Ixminus)=X(1,Ixminus)+2*dim;
    X(2,Iyplus)=X(2,Iyplus)-2*dim;
    X(2,Iyminus)=X(2,Iyminus)+2*dim;
    X(3,Izplus)=X(3,Izplus)-2*dim;
    X(3,Izminus)=X(3,Izminus)+2*dim;
    
    %Force and field calculations

    rx  = squareform(pdist(X',diffx));
    ry  = squareform(pdist(X',diffy));
    rz  = squareform(pdist(X',diffz));
    
    rx  = (rx.*corrector)';
    ry  = (ry.*corrector)';
    rz  = (rz.*corrector)';
  
    r   = (rx.^2+ry.^2+rz.^2).^(1/2)+eye(N);
    phi = r/2/dim;
    
    rhatx = rx./r;
    rhaty = ry./r;
    rhatz = rz./r;
    
    % Declare the time step
    %h = abs(t-i*0.0001);

%     Dnx = ones(N,1)*D(1,:);
%     Dny = ones(N,1)*D(2,:);
%     Dnz = ones(N,1)*D(3,:);
% 
%     Dmx = D(1,:)'*ones(1,N);
%     Dmy = D(2,:)'*ones(1,N);
%     Dmz = D(3,:)'*ones(1,N);
% 
%     DnR = rhatx.*Dnx+rhaty.*Dny+rhatz.*Dnz;
%     DmR = rhatx.*Dmx+rhaty.*Dmy+rhatz.*Dmz;
% 
%     DmDn = Dnx.*Dmx + Dny.*Dmy+ Dnz.*Dmz;

    %Fmagx = (3*Q)*((1./r.^4-eye(N)).*(rhatx.*DmDn-5*rhatx.*DmR.*DnR+Dnx.*DmR+Dmx.*DnR));
    %Fmagy = (3*Q)*((1./r.^4-eye(N)).*(rhaty.*DmDn-5*rhaty.*DmR.*DnR+Dny.*DmR+Dmy.*DnR));
    %Fmagz = (3*Q)*((1./r.^4-eye(N)).*(rhatz.*DmDn-5*rhatz.*DmR.*DnR+Dnz.*DmR+Dmz.*DnR));

    % ELECTRIC FIELD
    % In this section, the local field directions are found by determining
    % the local coordinates of each particle to write all local dipole
    % directions, and transform back into global coordinates

    % First use transformatins to determine local x,y,z coordinates
%     xold = [1,0,0];
%     yold = [0,1,0];
%     zold = [0,0,1];
% 
%     znew = D(:,:,1);
%     xnew = zeros(N,3);
%     ynew = zeros(N,3);
% 
%     rzq = zeros(N,4);
% 
%     for iz = 1:N
% 
%         rzq(iz,:) = vrrotvec(zold,D(:,iz,1)); % find rotation axis and angle
% 
%         xnew(iz,:) = rodrigues_rot(xold,rzq(iz,1:3),rzq(iz,4)); % find local x-axis
%         ynew(iz,:) = rodrigues_rot(yold,rzq(iz,1:3),rzq(iz,4)); % find local y-axis
% 
%     end

    % Next, dot product E0 with Dx, Dy, Dz for change in bases
%     E0nx = dot(xnew,E',2); % dotproduct of dx with E0  - Field in the local x-axis of the particle (coincides with a)
%     E0ny = dot(ynew,E',2); % Dot product of Dy with E0 -Field in the local y-axis of the particle (coincides with b)
%     E0nz = dot(znew,E)'; % Dot product of Dy with E0 - Field in the local z-axis of the particle (coincides with c)
% 
%     plocx = V*(eps2 - eps1).*E0nx./(1 + ((eps2-eps1)/eps1)*Ly);
%     plocy = V*(eps2 - eps1).*E0ny./(1 + ((eps2-eps1)/eps1)*Ly);
%     plocz = V*(eps2 - eps1).*E0nz./(1 + ((eps2-eps1)/eps1)*Ly);
% 
%     pglobx = plocx.*xnew(:,1) + plocy.*ynew(:,1) + plocz.*znew(1,:)';
%     pgloby = plocx.*xnew(:,2) + plocy.*ynew(:,2) + plocz.*znew(2,:)';
%     pglobz = plocx.*xnew(:,3) + plocy.*ynew(:,3) + plocz.*znew(3,:)';
% 
%     pmR = rhatx.*pglobx + rhaty.*pgloby + rhatz.*pglobz;
% 
%     Eddx = (1/4/pi/eps0)*sum((1./r.^3-eye(N)).*(3*pmR.*rhatx-pglobx)); % dipole-induced field effects in x-dir
%     Eddy = (1/4/pi/eps0)*sum((1./r.^3-eye(N)).*(3*pmR.*rhaty-pgloby)); % dipole-induced field effects in y-dir
%     Eddz = (1/4/pi/eps0)*sum((1./r.^3-eye(N)).*(3*pmR.*rhatz-pglobz)); % dipole-induced field effects in z-dir

    %pmx = pglobx*ones(1,N);
    %pmy = pgloby*ones(1,N); 
    %pmz = pglobz*ones(1,N);

    %dEdx = sum((rx.*(-3).*1/4/pi/eps0).*(1./r.^(5)-eye(N)).*(3*pmR.*rhatx-pmx));
    %dEdy = sum((ry.*(-3).*1/4/pi/eps0).*(1./r.^(5)-eye(N)).*(3*pmR.*rhaty-pmy));
    %dEdz = sum((rz.*(-3).*1/4/pi/eps0).*(1./r.^(5)-eye(N)).*(3*pmR.*rhatz-pmz));

    %E = [Ex; Ey; Ez];

    % DEP calculations - PROBLEM FOUND: THESE SHOULD BE WRITTEN IN LOCAL
    % COORDINATES, AND THEN GLOBAL
    %Fdepx = V.*(eps2 - eps1).*Ex.*dEdx./(1 + ((eps2-eps1)/eps1).*Lx); % DEP force in x direction
    %Fdepy = V.*(eps2 - eps1).*Ey.*dEdy./(1 + ((eps2-eps1)/eps1).*Ly); % DEP force in y direction
    %Fdepz = V.*(eps2 - eps1).*Ez.*dEdz./(1 + ((eps2-eps1)/eps1).*Lz); % DEP force in z direction

    %Fdepx = pmx.*dEdx;
    %Fdepy = pmx.*dEdy;
    %Fdepz = pmx.*dEdz;
    %pnx = ones(N,1)*pglobx';
    %pny = ones(N,1)*pgloby';
    %pnz = ones(N,1)*pglobz';

    %pmx = pglobx*ones(1,N);
    %pmy = pgloby*ones(1,N);
    %pmz = pglobz*ones(1,N);

    %pnR = rhatx.*pnx+rhaty.*pny+rhatz.*pnz;
    %pmR = rhatx.*pmx+rhaty.*pmy+rhatz.*pmz;

    %pmpn = pnx.*pmx + pny.*pmy+ pnz.*pmz;

    %Fdepx = (1/(4*pi*eps0)).*(3)*((1./r.^4-eye(N)).*(rhatx.*pmpn-5*rhatx.*pmR.*pnR+pnx.*pmR+pmx.*pnR));
    %Fdepy = (1/(4*pi*eps0)).*(3)*((1./r.^4-eye(N)).*(rhaty.*pmpn-5*rhaty.*pmR.*pnR+pny.*pmR+pmy.*pnR));
    %Fdepz = (1/(4*pi*eps0)).*(3)*((1./r.^4-eye(N)).*(rhatz.*pmpn-5*rhatz.*pmR.*pnR+pnz.*pmR+pmz.*pnR));

    % Repulsive forces
    %Frepx = (QQ)*rhatx.*exp(-beta*((r+corrector2)/2/R-1));
    %Frepy = (QQ)*rhaty.*exp(-beta*((r+corrector2)/2/R-1));
    %Frepz = (QQ)*rhatz.*exp(-beta*((r+corrector2)/2/R-1));

    %Fx = sum(Fmagx + Frepx + Fdepx);
    %Fy = sum(Fmagy + Frepy + Fdepy);
    %Fz = sum(Fmagz + Frepz + Fdepz);

    % Magnetic Field calculations
    %Hddx = (M*V/4/pi)*sum((1./r.^3-eye(N)).*(3*DmR.*rhatx-Dmx));
    %Hddy = (M*V/4/pi)*sum((1./r.^3-eye(N)).*(3*DmR.*rhaty-Dmy));
    %Hddz = (M*V/4/pi)*sum((1./r.^3-eye(N)).*(3*DmR.*rhatz-Dmz));

    %Hx = H0x*ones(1,N) + Hddx;
    %Hy = H0y*ones(1,N) + Hddy;
    %Hz = H0z*ones(1,N) + Hddz;

    % The dot product of D and H
    %DnH = Dnx(1,:).*Hx + Dny(2,:).*Hy+ Dnz(2,:).*Hz;

    %Eh = -mu0*M*V*(Hx.*D(1,:)+Hy.*D(2,:)+Hz.*D(3,:));

    % Calculate angular velocity due to magnetic field assuming zero angular acceleration
    %Wx = (QQQ)*(Hx-Dnx(1,:).*DnH);
    %Wy = (QQQ)*(Hy-Dny(1,:).*DnH);
    %Wz = (QQQ)*(Hz-Dnz(1,:).*DnH);

    % Electric torque calculations

    % First calculate field at each particle location
%     Ex = E0x*ones(1,N) + Eddx; % total field in x
%     Ey = E0y*ones(1,N) + Eddy; % total field in y
%     Ez = E0z*ones(1,N) + Eddz; % total field in z

%    Eloc = [Ex; Ey; Ez];

    % Transform field from global to local
%     Enx = dot(xnew,Eloc',2); % dotproduct of dx with E0  - Field in the local x-axis of the particle (coincides with a)
%     Eny = dot(ynew,Eloc',2); % Dot product of Dy with E0 -Field in the local y-axis of the particle (coincides with b)
%     Enz = dot(znew,Eloc)'; % Dot product of Dy with E0 - Field in the local z-axis of the particle (coincides with c)

    % Determine local torques for each particle
%     Tex = V*((eps2 - eps1)^2).*(Lz - Ly).*Eny.*Enz./(eps1*(1 + ((eps2-eps1)/eps1)*Ly)*(1 + ((eps2-eps1)/eps1)*Lz));
%     Tey = V*((eps2 - eps1)^2).*(Lx - Lz).*Enx.*Enz./(eps1*(1 + ((eps2-eps1)/eps1)*Lx)*(1 + ((eps2-eps1)/eps1)*Lz));
%     Tez = V*((eps2 - eps1)^2).*(Ly - Lx).*Enx.*Eny./(eps1*(1 + ((eps2-eps1)/eps1)*Lx)*(1 + ((eps2-eps1)/eps1)*Ly));

    % Global torques
%     Tnewx = Tex.*xnew(:,1) + Tey.*ynew(:,1) + Tez.*znew(1,:)';
%     Tnewy = Tex.*xnew(:,2) + Tey.*ynew(:,2) + Tez.*znew(2,:)';
%     Tnewz = Tex.*xnew(:,3) + Tey.*ynew(:,3) + Tez.*znew(3,:)';

    % Calculate T x d so that it can be added into the angular EoM
    %Wex = QQQe.*(Tnewy'.*Dnz(1,:) - Tnewz'.*Dny(1,:)); % ELECTRIC TORQUE IN X - DIR 
    %Wey = QQQe.*(Tnewz'.*Dnx(1,:) - Tnewx'.*Dnz(1,:)); % ELECTRIC TORQUE IN Y - DIR 
    %Wez = QQQe.*(Tnewx'.*Dny(1,:) - Tnewy'.*Dnx(1,:)); % ELECTRIC TORQUE IN Z - DIR 

    %Tx = Wx + Wex;
    %Ty = Wy + Wey;
    %Tz = Wz + Wez;

    % Set the magnetic field based on the current time
%     H0z = 10^5*sin((t+timeScale*h)*4*pi/10);

    % Magnetic Field calculations
%     Hddx = (M*V/4/pi)*sum((1./r.^3-eye(N)).*(3*DmR.*rhatx-Dmx));
%     Hddy = (M*V/4/pi)*sum((1./r.^3-eye(N)).*(3*DmR.*rhaty-Dmy));
%     Hddz = (M*V/4/pi)*sum((1./r.^3-eye(N)).*(3*DmR.*rhatz-Dmz));

%     Hx = H0x*ones(1,N) + Hddx;
%     Hy = H0y*ones(1,N) + Hddy;
%     Hz = H0z*ones(1,N) + Hddz;

    % Calculate magnetic torque on each particle
%     Tmx = mu0*M*V/J.*(Dny(1, :).*Hz - Dnz(1, :).*Hy);
%     Tmy = mu0*M*V/J.*(Dnz(1, :).*Hx - Dnx(1, :).*Hz);
%     Tmz = mu0*M*V/J.*(Dnx(1, :).*Hy - Dny(1, :).*Hx);

    % Calculate the drag torque on each particle
%     Tdx = 8*pi*eta*R^3/J.*omega(1, :);
%     Tdy = 8*pi*eta*R^3/J.*omega(2, :);
%     Tdz = 8*pi*eta*R^3/J.*omega(3, :);

    % Calculate the total torque on each particle
%     Tx = Tmx - Tdx + Tnewx;
%     Ty = Tmy - Tdy + Tnewy;
%     Tz = Tmz - Tdz + Tnewz;

    %vx = Fx/Ds;
    %vy = Fy/Ds;
    %vz = Fz/Ds;

    %v = [vx;vy;vz];
%    T = [Tx;Ty;Tz];

    %Numerical Integration

    %We determine the value of the time step:
    %h = min([10^-1/max(max(T)),abs(t-i*0.0001)]);

    % Correct repulsive force: create array of forces dependent on the
    % relative magnitude each particle is experiencing?


    %h = min([10^-7/max(max(v)),10^-5/max(max(T)),abs(t-i*0.0001)]);

    %We compute the evolution of the states:
    %X = X + v*h;%UPDATES POSITION
    %D = D + T*h;
    
    init = (horzcat(omega(1:3*N), D(1:3*N)))';
    
    [t, data] = ode45(@(t, params) odefcn(t, params, N, rhatx, rhaty, rhatz, r, E, E0x, E0y, E0z, M, V, eps0, eps1, eps2, mu0, Lx, Ly, Lz, H0x, H0y, J, eta, R), [0 0.001], init);
    
    [it, ~] = size(t);
    
    finalD = data(end, (3*N + 1):6*N);
    
    %Calculate the new D
    for i = 1:3*N
        if rem(i - 1, 3) == 0
            D(1, floor(i/3) + 1) = finalD(end, i);
        elseif rem(i - 2, 3) == 0
            D(2, floor(i/3) + 1) = finalD(end, i);
        else
            D(3, floor(i/3)) = finalD(end, i);
        end
    end
    
    %D = data((3*N + 1):6*N);
    D = D./(ones(3,1)*sum(D.^2).^(1/2));
    %t  = t+h;
    
    plotcond = 10;
    
    for i = 1:plotcond:it
        for j = 1:3*N
            if rem(j - 1, 3) == 0
                tempD(1, floor(j/3) + 1) = data(i, j);
            elseif rem(j - 2, 3) == 0
                tempD(2, floor(j/3) + 1) = data(i, j);
            else
                tempD(3, floor(j/3)) = data(i, j);
            end
        end
        
        [avg(floor(i/plotcond) + 1), std(floor(i/plotcond)+1)] = zAngles(tempD);
        t(i)
    end
    
    avg
    std
    
%     if rem(it, 500) == 0 && it <= 2000
%         Ddata(:, :, it/500 + 1) = D;
%     end
%    maxDeflection = 180/pi*max(acos(sum(D(:,:).*Ddata(:,:,1), 2)./(sqrt(sum(D(:,:).^2,2)).*sqrt(sum(Ddata(:,:,1).^2,2)))));
    
        %Increasing the loop parameters
%    it = it+1;

        % Keeping track of running time
%    tf = cputime-ti;

%     plotcond = rem(it,plotmult); 
%     if plotcond == 0
%         it
%        t
%     end
    
%     plotcond = rem(it,plotmult); 
%      
%     if plotcond == 0
%         
%         filename = sprintf('Dynsim_N=%d_Ex=%d_Ey=%d_Ez=%d_Hz=%d_maxit%d_it%d',N,E0x,E0y,E0z,H0z,maxit,it);       
%         save(filename,'X','D','times')
%         
%         zoldplot = [0,0,1];
% 
%         
%         ff = figure('visible','off');
%         hold on
%             for ipl = 1:N
% 
%                 [xplot, yplot, zplot] = ellipsoid(X(1,ipl),X(2,ipl),X(3,ipl),a,b,c,20);
% 
%                 S(ipl) = surf(xplot, yplot, zplot);%,'EdgeColor','none');
% 
%                 caxis('manual');
% 
%                 S(ipl).FaceLighting = 'gouraud';
% 
%                 S(ipl).EdgeColor = 'none';
%                 %colormap copper
%                 %caxis([50 80]*10^-7)
% 
%                 rplot(ipl,:) = vrrotvec(zoldplot,D(:,ipl));
% 
%                 rotate(S(ipl),rplot(ipl,1:3),rplot(ipl,4)*(180/pi),[X(1,ipl),X(2,ipl),X(3,ipl)])
% 
%             end
%             
%             %light('Position',[0 0 6*10^-6],'Style','local')
%             light()
%         
%             quiver3(X(1,:),X(2,:),X(3,:),D(1,:),D(2,:),D(3,:))
%             quiver3(X(1,:),X(2,:),X(3,:),pglobx',pgloby',pglobz')
%             view(130,40)
%                       
%             str3 = sprintf('N = %0.2d, t = %0.1ds, alpha = %0.1d,\ntime steps: %d, h = %d',N, t, vol, it,h);
%             title(str3)
% 
%             axis([-dim*1.02 dim*1.02 -dim*1.02 dim*1.02 -dim*1.02 dim*1.02])
%             xlabel('[m]')
%             ylabel('[m]')
%             zlabel('[m]')
%             filename = sprintf('Dynsim_N=%d_Ex=%d_Ey=%d_Ez=%d_Hz=%d_maxit%d_it%d.png',N,E0x,E0y,E0z,H0z,maxit,it);       
%  
%             saveas(ff,filename)
%             
%         
%         %maxFdep = max([max(Fdepx), max(Fdepy), max(Fdepz)])
%         %maxFmag = max([max(sum(Fmagx)), max(sum(Fmagy)), max(sum(Fmagz))])
%         
%         %maxWe = max([max(Wex), max(Wey), max(Wez)])
%         %maxWm = max([max(Wx), max(Wy), max(Wz)])
%         
%         toc
%     end
    
%end

%     tdata(i) = t;
%     
%     DDdata(:,:) =  ((1./phi.^3-eye(N)).*(3*DmR.*rhatx-Dmx)).*Dnx...
%                     +((1./phi.^3-eye(N)).*(3*DmR.*rhaty-Dmy)).*Dny...
%                     +((1./phi.^3-eye(N)).*(3*DmR.*rhatz-Dmz)).*Dnz;
%     
%     f(i) = (1/4/pi/N^2)*sum(sum(DDdata(:,:)));
%     g(i) = (1/N)*sum(D(3,:));
%     
%     i = i+1;  
%end
    Xdata(:,:,2)  = X;
    Ddata(:,:,2)  = D;
    
    %maxDeflection1 = 180/pi*max(acos(sum(Ddata(:,1:2,2).*Ddata(:,1:2,1), 2)./(sqrt(sum(Ddata(:,1:2,2).^2,2)).*sqrt(sum(Ddata(:,1:2,1).^2,2)))))
    %maxDeflection2 = 180/pi*max(acos(sum(Ddata(:,2:3,2).*Ddata(:,2:3,1), 2)./(sqrt(sum(Ddata(:,2:3,2).^2,2)).*sqrt(sum(Ddata(:,2:3,1).^2,2)))))
    %maxDeflection = 180/pi*max(acos(sum(Ddata(:,:,2).*Ddata(:,:,1), 2)./(sqrt(sum(Ddata(:,:,2).^2,2)).*sqrt(sum(Ddata(:,:,1).^2,2)))))
    
    H0 = [H0x; H0y; H0z];
    E0 = [E0x; E0y; E0z];

    %maxFmag = [max(sum(Fmagx)), max(sum(Fmagy)), max(sum(Fmagz))]
    %maxFrep = [max(sum(Frepx)), max(sum(Frepy)), max(sum(Frepz))];
    %maxFdep = [max(sum(Fdepx)), max(sum(Fdepy)), max(sum(Fdepz))]
    
%     size(Wx);
%     maxTmag = [max(Wx), max(Wy), max(Wz)]
%     size(Wex)
%     maxTe = [max(Wex), max(Wey), max(Wez)]
   
data=struct('Parameters',[N,vol,HtoM,R,dim,M,eta],'X',Xdata,'D',Ddata,'g',g,'f',f,'t1',tdata,'E0',E0,'H0',H0,'dim',dim);
%,'pglobx',pglobx,'pgloby',pgloby,'pglobz',pglobz);

end