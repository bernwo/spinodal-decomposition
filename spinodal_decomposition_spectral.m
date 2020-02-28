%% CLEAR EVERYTHING
clear
close all
clc
magma_cMap = getPyPlot_cMap('magma');
%% Initialisation
close all
steps               =    1000;
save_interval       =    50;
mu                  =    0; %mean
sigma               =    1e-3; %standard deviation
gamma               =    0.5;
nx                  =    128*4; % number of cells in the x-direction
ny                  =    128*4; % number of cells in the y-direction
dx                  =    1.0;
dy                  =    1.0;
dt                  =    0.8;
k2                  =    zeros(ny,nx);
Lx                  =    nx*dx;
Ly                  =    ny*dy;

% Initialising k^2 space
for i = 1:nx
    for j = 1:ny
        if ((i-1) < nx/2)
           kx2     =    (2*pi*(i)/Lx)^2;
        else
           kx2     =    (2*pi*(i-nx)/Lx)^2;
        end
        if ((j-1) < ny/2)
           ky2     =    (2*pi*(j)/Ly)^2;
        else
           ky2     =    (2*pi*(j-ny)/Ly)^2;
        end
        
        k2(j,i)     =    kx2 + ky2;
    end
end
% Visualising the k^2 space
% imagesc(k2);
% axis equal
% xlim([0 nx]);
% ylim([0 ny]);
% colorbar()

%% Animating
filename = 'spinodal.gif';
c                   =    normrnd(mu,sigma,[ny,nx]);
firsttime           =    true;

% Begin timestep
for nn = 1:steps
    ck                  =    fft2(c);
    nonlinear           =    fft2(-c.^3 + c);
    ck                  =    (ck + k2.*dt.*nonlinear) ./ (k2.^2.*gamma.*dt+1);
    c                   =    real(ifft2(ck));
    % Write to the GIF File 
    if(mod(nn,save_interval) == 0)
        disp(nn)
        h = figure('Renderer', 'painters', 'Position', [1 1 700 700]);
        imagesc(real(c));
        axis equal
        axis off tight
        xlim([0 nx]);
        ylim([0 ny]);
%         colorbar()
        colormap(magma_cMap)
        frame = getframe(gca); 
%         im = frame2im(frame); 
%         [imind,cm] = rgb2ind(im,256);
        [imind,cm] =  rgb2ind(frame.cdata,256,'nodither');
        if firsttime
            imwrite(imind,cm,filename,'gif','DelayTime',0.05,'Loopcount',inf);
            firsttime       =    false;
        else
            imwrite(imind,cm,filename,'gif','DelayTime',0.05,'WriteMode','append'); 
        end
        close
    end
end