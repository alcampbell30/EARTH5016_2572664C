%*****  1D ADVECTION DIFFUSION MODEL OF HEAT TRANSPORT  *******************

%*****  Initialise Model Setup

% create x-coordinate vectors
xc = h/2:h:W-h/2;      % x-coordinate vector for cell centre positions [m]
zc = h/2:h:D-h/2;      % z-coordinate vector for cell centre positions [m]
xf = 0:h:W;            % x-coordinate vectore for cell face positions [m]
zf = 0:h:D;            % z-coordinate vectore for cell face positions [m]

[Xc,Zc] = meshgrid(xc,zc);  % create 2D coordinate arrays

% set up index array for boundary conditions
ix3 = [          Nx,1:Nx,1       ];  % periodic sides
ix5 = [     Nx-1,Nx,1:Nx,1 ,2    ];
ix7 = [Nx-2,Nx-1,Nx,1:Nx,1 ,2 ,3 ];

iz3 = [           1,1:Nz,Nz      ];  % closed/insulating top/bot
iz5 = [        1, 1,1:Nz,Nz,Nz   ];
iz7 = [   1,   1, 1,1:Nz,Nz,Nz,Nz];

% create smooth random perturbation field
rng(15);
dr  = randn(Nz,Nx);

for ii = 1:10
    dr  = dr + (diff(dr(iz3,:),2,1) + diff(dr(:,ix3),2,2))/8;
end

% set initial condition for temperature at cell centres
T   = Ttop + (Tbot-Ttop)./D.*Zc + dr*1;  % initialise T array on linear gradient

% initialise density and mobility
rho = rho0.*(1 - aT.*(T-Ttop));
KD  = KD0.*ones(Nz,Nx);
kT  = kT0.*ones(Nz,Nx);

% initialise Darcy velocity and non-hydrostatic pressure
u = zeros(Nz,Nx+1);
w = zeros(Nz+1,Nx);
p = zeros(Nz,Nx);

% initialise residual and update fields
res_T = 0*T;  upd_T = 0*T;  dTdt = 0*T;
res_p = 0*p;  upd_p = 0*p;

% initialise output figure with initial condition
figure(1); clf
makefig(xc,zc,T,p,w,u,0,yr)

%*****  Solve Model Equations

dt = CFL * min((h/2)/max([w(:);u(:)]),(h/2)^2/max(kT(:))); % initial time step [s]
t  = 0;  % initial time [s]
k  = 0;  % initial time step count

% loop through time steps until stopping time reached
while t <= tend

    % increment time and step count
    t = t+dt;
    k = k+1;

    % print time step header
    fprintf(1,'\n\n*****  step = %d;  dt = %1.3e;  time = %1.3e \n\n',k,dt,t)

    % store old temperature and rate
    dTdto = dTdt;
    To    = T;

    % get time step step size
    dt    = CFL * min((h/2)/max([w(:);u(:)]),(h/2)^2/max(kT(:))); % time step [s]

    resnorm = 1;  % initialise residual norm
    it      = 0;  % initialise iteration count

    % loop through pseudo-transient iterations until convergence criterion reached
    while resnorm > tol

        % update temperature every 'nup' iterations
        if ~mod(it,nup) && k>1

            % get T-dependent segregation mobility
            kT  = kT0 + cT.*max(0,T).^mT;

            % get rate of change
            dTdt = diffusion(T,kT,h,ix3,iz3) + advection(T,u,w,h,ix7,iz7,ADVN);

            % get temperature residual
            res_T = (T - To)/dt - (dTdt + dTdto)/2;

            % set isothermal boundaries on top/bot
            res_T(1  ,:) = 0;
            res_T(end,:) = 0;

            % get solution update
            upd_T = - alpha*res_T*dt/2;

            % update solution
            T     = T + upd_T;

        end

        % get density and density contrast
        rho   = rho0.*(1 - aT.*(T-Ttop));  % T-dependent density
        Drho  = rho - mean(rho,2);         % subtract horizontal mean
        Drhoz = (Drho(iz3(1:end-1),:)+Drho(iz3(2:end),:));  % on z-faces
        Drhoz([1 end],:) = 0;              % no flow across top/bot bounds

        % get p-dependent segregation mobility
        KD  = KD0 + cp.*max(0,p).^mp;
        KDx = (complete);  % on x-faces
        KDz = (complete);  % on z-faces

        % get Darcy flux (vD = - KD (Grad(p) - D(rho) g)
        u   = (complete);  % x-speed
        w   = (complete);  % z-speed

        % get pressure residual (Div.v = 0)
        res_p = (complete);

        % get pseudo-time step size
        dtau  = (complete);

        % get solution update
        upd_p = (complete);

        % update solution
        p     = (complete);

        it = it+1; % increment iteration count

        % get residual norm and print convergence every 'nup' iterations
        if ~mod(it,nup)
            resnorm = norm(upd_T(:),2)./norm(T(:)+eps,2) ...
                    + norm(upd_p(:),2)./norm(p(:)+eps,2);
            if isnan(resnorm); error('!!! Solver failed with nan !!!'); end
            fprintf(1,'     it = %d;  res = %e \n',it,resnorm); 
        end

    end

    % plot model progress every 'nop' time steps
    if ~mod(k,nop)
        makefig(xc,zc,T,p,w,u,t,yr);
    end

end


%*****  Utility Functions  ************************************************

% Function to make output figure
function makefig(x,z,T,p,w,u,t,yr)

clf; 

% plot temperature in subplot 1
subplot(2,2,1);
imagesc(x,z,T); axis equal tight; colorbar; hold on
contour(x,z,T,[100,150,200],'k');

ylabel('z [m]','FontSize',15)
title('Temperature [C]','FontSize',17)

% plot pressure in subplot 2
subplot(2,2,2)
imagesc(x,z,p); axis equal tight; colorbar

xlabel('x [m]','FontSize',15)
ylabel('z [m]','FontSize',15)
title('Pressure [Pa]','FontSize',17)

% plot z-speed in subplot 3
subplot(2,2,3)
imagesc(x,z,-w*yr); axis equal tight; colorbar

xlabel('x [m]','FontSize',15)
ylabel('z [m]','FontSize',15)
title('z-Speed [m/yr]','FontSize',17)

% plot x-speed in subplot 1
subplot(2,2,4)
imagesc(x,z,u*yr); axis equal tight; colorbar

xlabel('x [m]','FontSize',15)
ylabel('z [m]','FontSize',15)
title('x-Speed [m/yr]','FontSize',17)

sgtitle(['time = ',num2str(t/yr),' [yr]'],'FontSize',17)
drawnow;

end

% Function to calculate diffusion rate
function [dTdt] = diffusion(f,k,h,ix,iz)
% calculate heat flux by diffusion

qx = - k .* diff(f(ix, :), 1, 2)/h;
qz = - k .* diff(f(iz, :), 1, 1)/h;

% calculate flux balance for rate of change
dTdt = - (diff(qx)/h+diff(qz)/h);


end


% Function to calculate advection rate
function dTdt = advection(f,u,w,h,ix,iz,ADVN)

% split the velocities into positive and negative
upos = 0.5*(u + abs(u));    % positive velocity in x direction
uneg = 0.5*(u - abs(u));    % negative velocity in x direction

wxpos = 0.5*(w + abs(w)); %positive velocity in the z direction
wneg = 0.5*(w-abs(w)); %negative velocity in the z direction

% get values on stencil nodes for x direction and velocity 
fxmmm = f(ix(1:end-6)); %variable (temperature in this case) value in the cell 3 to the left of central cell (m = minus)
fxmm  = f(ix(2:end-5)); %variable (temperature in this case) value in the cell 2 to the left of central cell
fxm   = f(ix(3:end-4)); %variable (temperature in this case) value in the cell 1 to the left of central cell
fxc   = f(ix(4:end-3)); %variable (temperature in this case) value in the central cell (c = central)
fxp   = f(ix(5:end-2)); 
fxpp  = f(ix(6:end-1));
fxppp = f(ix(7:end-0));

fzmmm = f(iz(1:end-6)); %variable (temperature in this case) value in the cell 3 to the left of central cell (m = minus)
fzmm  = f(iz(2:end-5)); %variable (temperature in this case) value in the cell 2 to the left of central cell
fzm   = f(iz(3:end-4)); %variable (temperature in this case) value in the cell 1 to the left of central cell
fzc   = f(iz(4:end-3)); %variable (temperature in this case) value in the central cell (c = central)
fzp   = f(iz(5:end-2)); %variable (temperature in this case) value in the cell 1 to the right of the central cell (p = plus)
fzpp  = f(iz(6:end-1));
fzppp = f(iz(7:end-0));

% calculate heat flux by advection
switch ADVN
    case 'UPW1'
        fxppos = fxc;     fxpneg = fxp;
        fxmpos = fxm;     fxmneg = fxc;

        fzppos = fzc;     fzpneg = fzp;
        fzmpos = fzm;     fzmneg = fzc;

    case 'CFD2'
        fxppos = (fxc+fxp)./2;     fxpneg = fxppos;
        fxmpos = (fxc+fxm)./2;     fxmneg = fxmpos;

        fzppos = (fzc+fzp)./2;   fzpneg = fzp;
        fzmpos = (fzc+fzm)./2;     fzmneg = fzc;

    case 'FRM2'
        fxppos = fxc + (fxp-fxm )./4;     fxpneg = fxp + (fxc-fxpp)./4;
        fxmpos = fxm + (fxc-fxmm)./4;     fxmneg = fxc + (fxm-fxp )./4;
    
        fzppos = fzc + (fzp-fzm )./4;     fzpneg = fzp + (fzc-fzpp)./4;
        fzmpos = fzm + (fzc-fzmm)./4;     fzmneg = fzc + (fzm-fzp )./4;

    case 'UPW3'
        fxppos = (2*fxp + 5*fxc - fxm )./6;     fxpneg = (2*fxc + 5*fxp - fxpp)./6;
        fxmpos = (2*fxc + 5*fxm - fxmm)./6;     fxmneg = (2*fxm + 5*fxc - fxp )./6;

        fzppos = (2*fzp + 5*fzc - fzm )./6;     fzpneg = (2*fzc + 5*fzp - fzpp)./6;
        fzmpos = (2*fzc + 5*fzm - fzmm)./6;     fzmneg = (2*fzm + 5*fzc - fzp )./6;

    case 'WENO5'
        fxppos = weno5poly(fxmm ,  fxm, fxc, fxp, fxpp);     fxpneg = weno5poly(fxppp, fxpp, fxp, fxc, fxm );
        fxmpos = weno5poly(fxmmm, fxmm, fxm, fxc, fxp );     fxmneg = weno5poly( fxpp,  fxp, fxc, fxm, fxmm);

        fzppos = weno5poly(fzmm ,  fzm, fzc, fzp, fzpp);     fzpneg = weno5poly(fzppp, fzpp,fzp, fzc, fzm );
        fzmpos = weno5poly(fzmmm, fzmm, fzm, fzc, fzp );     fzmneg = weno5poly( fzpp, fzp, fzc, fzm, fzmm);
end

% calculate flux balance for rate of change in x direction
div_qxpos = upos.*(fxppos - fxmpos)/h;
div_qxneg = uneg.*(fxpneg - fxmneg)/h;
div_qx    = div_qxpos + div_qxneg;

div_qzpos = wpos.*(fzppos - fzmpos)/h;
div_qzneg = wneg.*(fzpneg - fzmneg)/h;
div_qz    = div_qzpos+div_qzneg;

dTdt = - (div_qx+div_qz);


end


% Function to calculate 5th order WENO polynomials from Jiang & Shu, 1996, J Comp Physics
function [fface] = weno5poly (fmm, fm, fc, fp, fpp)

% 5th order polynomials
p1 = (2*fmm - 7*fm + 11*fc )/6;
p2 = ( -fm  + 5*fc +  2*fp )/6;
p3 = (2*fc  + 5*fp -    fpp)/6;

% smoothness measure
b1 = 13/12*(fmm - 2*fm + fc ).^2 + 1/4*(  fmm - 4*fm + 3*fc ).^2;
b2 = 13/12*(fm  - 2*fc + fp ).^2 + 1/4*(  fm  -          fp ).^2;
b3 = 13/12*(fc  - 2*fp + fpp).^2 + 1/4*(3*fc  - 4*fp +   fpp).^2;

% weights
g   = [1/10, 6/10, 3/10];
wp1 = g(1)./(b1.^2 + eps);
wp2 = g(2)./(b2.^2 + eps);
wp3 = g(3)./(b3.^2 + eps);

% get flux (normalize weights at the same time)
fface = (wp1.*p1 + wp2.*p2 + wp3.*p3) ./ (wp1 + wp2 + wp3) ;

end