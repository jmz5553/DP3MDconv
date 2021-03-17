function DPM3Dconv_run_file_New_NoTranslation(defaultN,volFrac,saveFile, runs)
% Multi-field processing - dynamic simulations run script
% dim for N = 100: 8.6823e-06 --> 86 microns, one particle is 0.5 microns
% dim for N = 200: 1.0939e-05
% dim for N = 500: 1.4847e-05

tic

Scell = cell(4,1);

parfor i = 1:runs

Scell{i,1} = DPM3DconvE_NoTranslation_ode45(defaultN, volFrac)
toc

end

formatString = '%d ';
for i = 1:defaultN
   formatString = strcat(formatString,' %d ');
end

data = Scell{1,1};
 X = data.X;
 
resPos = strcat(saveFile,'_resPos.txt');
fid = fopen(resPos,'wt');
for ii = 1:size(X,1)
    fprintf(fid,'%g\t',X(ii,:));
    fprintf(fid,'\n');
end
fclose(fid);
 

 D = data.D;

resDir = strcat(saveFile,'_resDir.txt');
fid = fopen(resDir,'wt');
for ii = 1:size(D,1)
    fprintf(fid,'%g\t',D(ii,:));
    fprintf(fid,'\n');
end
fclose(fid);


 t1 = data.t1;
 f = data.f;
 g = data.g;
 %pglobx = data.pglobx;
 %pgloby = data.pgloby;
 %pglobz = data.pglobz;
 N = data.Parameters(1);

% For regular results, N = 100
%save('No_field_t0p75_N100')

% For convergence study
save('E_only_CASE2_t0p5_N500_1to4');

%%
% PLOT ELLIPSOIDS FROM CURRENT WORKSPACE DATA

figure('visible','on')
hold on
%conf = 1; % 1 for before, 2 for after
%f1 = scatter3(X(1,:,conf),X(2,:,conf),X(3,:,conf),250,'filled')
%f1.MarkerFaceAlpha = 1;
%quiver3(X(1,:,conf),X(2,:,conf),X(3,:,conf),D(1,:,conf),D(2,:,conf),D(3,:,conf))

for conf = 1:2
    figure(conf)
    hold on
    
    Xp = X(:,:,1);
    Dp = D(:,:,conf);

    a = 0.5*10^-6;
    b = 0.5*10^-6;
    c = 0.15*10^-6;

    R = 0.5*10^-6;
    V    = 4/3*pi*R^3; %particle volume [m^3]
    syms a
    eq1 = V == 4/3*pi*a^2*(1/3.5)*a;
    a = max(real(double(solve(eq1))));
    b = a;
    c = a/3.5;

    zold = [0,0,1];

    for i = 1:N

        [x, y, z] = ellipsoid(Xp(1,i),Xp(2,i),Xp(3,i),a,b,c,20);

        S(i) = surf(x, y, z);%,'EdgeColor','none');

        caxis('manual');

        S(i).FaceLighting = 'gouraud';

        S(i).EdgeColor = 'none';
        %colormap copper
        %caxis([50 80]*10^-7)

        r(i,:) = vrrotvec(zold,Dp(:,i));

        rotate(S(i),r(i,1:3),r(i,4)*(180/pi),[Xp(1,i),Xp(2,i),Xp(3,i)])

    end
    
    light()
end
%light('Position',[0 0 6*10^-6],'Style','local')
%light()

%f2 = scatter3(X(1,:,conf),X(2,:,conf),X(3,:,conf),250,'filled')
%f2.MarkerFaceAlpha = 1;

vol = volFrac; % volumetric density
realtime = cputime; % 
it = 40000

% str3 = sprintf('N = %0.2d, t = %0.1ds, alpha = %0.1d,\nime steps: %d, h = 5.446809e-06',N, t1(end), vol, it);
% title(str3)

% f = figure('visible','off');
% quiver3(Xp(1,:),Xp(2,:),Xp(3,:),Dp(1,:),Dp(2,:),Dp(3,:));
% view(130,40);
% savefig(f,saveFile);

%quiver3(Xp(1,:),Xp(2,:),Xp(3,:),pglobx',pgloby',pglobz')
%%

end