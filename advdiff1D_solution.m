%*****  1D ADVECTION DIFFUSION MODEL OF HEAT TRANSPORT  *******************

%*****  Initialise Model Setup

% create x-coordinate vectors
xc = dx/2:dx:W-dx/2;    % coordinate vector for cell centre positions [m]
xf = 0:dx:W;            % coordinate vectore for cell face positions [m]

% set time step size
dt = CFL * min((dx/2)/u0,(dx/2)^2/k0); % time step [s]

% set up index array for boundary conditions
switch BC
    case 'periodic'
        % example periodic indexing for N=4
        %   [4,1,2,3,4,1]    % 3-point stencil
        % [3,4,1,2,3,4,1,2]  % 5-point stencil
        ind3 = [        N,1:N,1    ];
        ind5 = [    N-1,N,1:N,1,2  ];
        ind7 = [N-2,N-1,N,1:N,1,2,3];
    case 'insulating'
        % example non-periodic indexing for N=4
        %   [1,1,2,3,4,4]      % 3-point stencil
        % [1,1,1,2,3,4,4,4]    % 5-point stencil
        ind3 = [    1,1:N,N    ];
        ind5 = [  1,1,1:N,N,N  ];
        ind7 = [1,1,1,1:N,N,N,N]; %having 1 and N facing each other makes the flux 0 across the boundary.
end

% set initial condition for temperature at cell centres
T   = T0 + dT .* exp(-(xc-W/2).^2/(4*wT^2));     % initialise T array at Tr, .* is element wise vector multiplication, .^ is element wise power
Tin = T;                                         % store initial condition
Ta  = T;                                         % initialise analytical solution

% initialise output figure with initial condition
figure(1); clf
makefig(xc,T,Tin,Ta,0)

%*****  Solve Model Equations

t = 0;  % initial time [s]
k = 0;  % initial time step count

% assemble coefficient matrix for time-dependent coefficients
At = speye(N,N)/dt; % sparse matrix with 1/dt on diagonal
% assemble coefficient matrix for space-dependent coefficients
j = []; i = []; a = []; % initialise index and value lists

j = [j,ind3(2:end-1)]; % equation numbers for internal nodes i=1:N
i = [i,ind3(1:end-2)]; % degree-of-freedom numbers for nodes i-1
a = [a,(-u0/2/dx - k0/dx^2).*ones(1,N)]; % coeff. values for nodes i-1


j = [j,ind3(2:end-1)]; % equation numbers for internal nodes i=1:N
i = [i,ind3(2:end-1)]; % degree-of-freedom numbers for nodes i
a = [a,(    +2k0/dx^2).*ones(1,N)]; % coeff. values for nodes i, 1/dt factor is taken care of by At

j = [j,ind3(2:end-1)]; % equation numbers for internal nodes i=1:N
i = [i,ind3(3:end)]; % degree-of-freedom numbers for nodes i+1
a = [a,(u0/2/dx -k0/dx^2.*ones(1,N)]; % coeff. values for nodes i+1 

% place values a at positions (j,i) in NxN sparse matrix Ax
Ax = sparse(j,i,a,N,N);

while t <= tend %while the time is still below the set end time 

    % increment time and step count
    t = t+dt;
    k = k+1;

    switch TINT
        case 'FE1'  % 1st-order Forward Euler time integration scheme
            
            dTdt = diffusion(T,k0,dx,ind3) + advection(T,u0,dx,ind7,ADVN);

            T = T + dTdt * dt;

        case 'RK2'  % 2nd-order Runge-Kutta time integration scheme
            
            dTdt1 = diffusion(T           ,k0,dx,ind3) + advection(T           ,u0,dx,ind7,ADVN);
            dTdt2 = diffusion(T+dTdt1*dt/2,k0,dx,ind3) + advection(T+dTdt1*dt/2,u0,dx,ind7,ADVN);

            T = T + dTdt2 * dt;

        case 'HE2'  % 2nd-order Heun's time integration scheme
            
            dTdt1 = diffusion(T         ,k0,dx,ind3) + advection(T         ,u0,dx,ind7,ADVN);
            dTdt2 = diffusion(T+dTdt1*dt,k0,dx,ind3) + advection(T+dTdt1*dt,u0,dx,ind7,ADVN);

            T = T + (dTdt1 + dTdt2)/2 * dt;

        case 'RK4'  % 4th-order Runge-Kutta time integration scheme
            
            dTdt1 = diffusion(T           ,k0,dx,ind3) + advection(T           ,u0,dx,ind7,ADVN);
            dTdt2 = diffusion(T+dTdt1/2*dt,k0,dx,ind3) + advection(T+dTdt1/2*dt,u0,dx,ind7,ADVN);
            dTdt3 = diffusion(T+dTdt2/2*dt,k0,dx,ind3) + advection(T+dTdt2/2*dt,u0,dx,ind7,ADVN);
            dTdt4 = diffusion(T+dTdt3  *dt,k0,dx,ind3) + advection(T+dTdt3  *dt,u0,dx,ind7,ADVN);

            T = T + (dTdt1 + 2*dTdt2 + 2*dTdt3 + dTdt4)/6 * dt;

    end

    % get analytical solution at time t
    wTt = sqrt(wT^2 + 1*k0*t);
    Ta  = T0 + dT .* wT./wTt .* exp(-(xc-W/2-u0*t).^2/(4*wTt^2)) + dT .* wT./wTt .* exp(-(xc+W/2-u0*t).^2/(4*wTt^2));

    % plot model progress
    if ~mod(k,nop)
        makefig(xc,T,Tin,Ta,t/yr);
    end

end

%*****  calculate numerical error norm
Err = norm(T - Ta,2)./norm(Ta,2);
disp(' ');
disp(['Advection scheme: ',ADVN]);
disp(['Time integration scheme: ',TINT]);
disp(['Numerical error = ',num2str(Err)]);
disp(' ');

%*****  Utility Functions  ************************************************

% Function to make output figure
function makefig(x,T,Tin,Ta,t)

plot(x,Tin,'k:',x,T,'r-',x,Ta,'k--','LineWidth',1.5); axis tight; box on;

xlabel('x [m]','FontSize',15)
ylabel('T [C]','FontSize',15)
title(['Temperature; time = ',num2str(t),' yr'],'FontSize',18)

drawnow;

end

% Function to calculate diffusion rate
function dTdt = diffusion(f,k0,dx,ind)

% calculate heat flux by diffusion
q = - k0 .* diff(f(ind))/dx;

% calculate flux balance for rate of change
dTdt = - diff(q)/dx;

end


% Function to calculate advection rate
function dTdt = advection(f,u0,dx,ind,ADVN)

% split the velocities into positive and negative
upos = 0.5*(u0 + abs(u0));    % positive velocity
uneg = 0.5*(u0 - abs(u0));    % negative velocity

% get values on stencil nodes
fmmm = f(ind(1:end-6));
fmm  = f(ind(2:end-5));
fm   = f(ind(3:end-4)); 
fc   = f(ind(4:end-3)); 
fp   = f(ind(5:end-2)); 
fpp  = f(ind(6:end-1));
fppp = f(ind(7:end-0));


% calculate heat flux by advection
switch ADVN
    case 'UPW1'
        fppos = fc;     fpneg = fp;
        fmpos = fm;     fmneg = fc;

    case 'CFD2'
        fppos = (fc+fp)./2;     fpneg = fppos;
        fmpos = (fc+fm)./2;     fmneg = fmpos;

    case 'FRM2'
        fppos = fc + (fp-fm )./4;     fpneg = fp + (fc-fpp)./4;
        fmpos = fm + (fc-fmm)./4;     fmneg = fc + (fm-fp )./4;

    case 'UPW3'
        fppos = (2*fp + 5*fc - fm )./6;     fpneg = (2*fc + 5*fp - fpp)./6;
        fmpos = (2*fc + 5*fm - fmm)./6;     fmneg = (2*fm + 5*fc - fp )./6;

    case 'WENO5'
        fppos = weno5poly(fmm ,  fm, fc, fp, fpp);     fpneg = weno5poly(fppp, fpp, fp, fc, fm );
        fmpos = weno5poly(fmmm, fmm, fm, fc, fp );     fmneg = weno5poly( fpp,  fp, fc, fm, fmm);
end

% calculate flux balance for rate of change
div_qpos = upos.*(fppos - fmpos)/dx;
div_qneg = uneg.*(fpneg - fmneg)/dx;
div_q    = div_qpos + div_qneg;

dTdt = - div_q;

end


function [fface] = weno5poly (fmm, fm, fc, fp, fpp)
% 5th order WENO polynomials from Jiang & Shu, 1996, J Comp Physics

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
