function [phi, r] = solver_CahnHilliard_2D_SAV_BDF(pde,domain,Nx,Ny,time,option)
% Cahn-Hilliard equation
% Reference:
%
%
% Qi Li
% Created: September 25, 2025

global dt kx ky kxx kyy k2 k4 hx hy Lx Ly ...
    epsilon M C0 S

if ~exist('option','var'), option = []; end

if ~isfield(option,'plotflag')
    option.plotflag = 0;
end
if ~isfield(option,'printflag')
    option.printflag = 0;
end
if ~isfield(option,'vtkflag')
    option.vtkflag = 0;
end
if ~isfield(option,'saveflag')
    option.saveflag = 0;
end
if ~isfield(option,'savefinal')
    option.savefinal = 0;
end
if ~isfield(option,'energyflag')
    option.energyflag = 0;
end
if ~isfield(option,'tol')
    option.tol = 10^-14;   % default tol
end
if ~isfield(option,'tolit')
    option.tolit = 10^-8;   % default tolit
end
if ~isfield(option,'maxit')
    option.maxit = 2000;   % default maxit
end

%%
T  = time.T;
t0  = time.t0;
dt = time.dt;
tsave = time.tsave;

dir_fig  = [pde.name '/fig'];
dir_data = [pde.name '/data'];

epsilon = pde.epsilon;
M       = pde.M;
C0      = pde.C0;
S       = pde.S;

Lx = domain.right - domain.left;
Ly = domain.top   - domain.bottom;
hx = Lx/Nx;
hy = Ly/Ny;
% x  = domain.left   + hx*(1:Nx);
% y  = domain.bottom + hy*(1:Ny);
x_c  = domain.left   + hx*(0:Nx-1);
y_c  = domain.bottom + hy*(0:Ny-1);
[xx, yy] = ndgrid(x_c,y_c);

% Fourier spectral 
% [k_x,k_y,kx,ky,kxx,kyy,k2,k4] = prepare_fft2_v2(Lx,Ly,Nx,Ny);
k_x = 1i*[0:Nx/2 -Nx/2+1:-1]*(2*pi/Lx);
k_y = 1i*[0:Ny/2 -Ny/2+1:-1]*(2*pi/Ly);
% [kx, ky] = ndgrid(k_x,k_y);
k2x = k_x.^2;
k2y = k_y.^2;
[kxx, kyy] = ndgrid(k2x,k2y);
k2 = kxx + kyy;
k4 = k2.^2;
% Furthermore, it is important to make the highest frequency N/2 to zero 
% in odd derivatives due to the symmetry.
k_x1 = 1i*[0:Nx/2-1 0 -Nx/2+1:-1]*(2*pi/Lx);
k_y1 = 1i*[0:Ny/2-1 0 -Ny/2+1:-1]*(2*pi/Ly);
[kx,  ky] = ndgrid(k_x1,k_y1);

phi_old = pde.initphi(xx,yy);

r_old = fun_r_init(phi_old);

nfigure =1;
nplot = round((T-t0)/dt);
nsave = round(tsave/dt);

tstart = tic;

option1 =  option;
option1.savefinal = 0;
time1 = time; time1.T = dt;
[phi_prev,r_prev] = solver_CahnHilliard_2D_SAV_1st(pde,domain,Nx,Ny,time1,option1);

for nt = 2:nplot
    t = t0 + nt*dt;

    phi_star = 2*phi_prev - phi_old;
    
    H = fun_H(phi_star);
    
    % Step 1 
    if isfield(pde,'exactphi') && isfield(pde,'rhsphi')
        tmp = pde.rhsphi(xx,yy,t) / M;
        % fprintf('aaa\n');
    else
        tmp = 0;
    end

    g1 = (4*r_prev - r_old)/3 - 1/2 * inner_int(H,(4*phi_prev-phi_old)/3);
    rhs = (4*phi_prev - phi_old) / (2*M*dt) + diff_lap(H).*g1 + tmp;

    % Step 2
    psi_B = inv_A(diff_lap(H));
    psi_C = inv_A(rhs);

    gamma =  -1/2 * inner_int(H,psi_B);
    Hphi = inner_int(H,psi_C) / (1 + gamma);

    % Step 3
    phi = 1/2 * Hphi .* psi_B + psi_C;

    r = 1/2*Hphi + g1;

    % update
    phi_old = phi_prev;  phi_prev = phi;
    r_old = r_prev;  r_prev = r;

    if 1 == option.energyflag
        calculate_energy1(out1,out2,hx,hy,t,phi10, phi20, phi30,u0);
    end

    if  0 == mod(nt,nsave)
        if 1 == option.printflag
            timeElapsed = toc(tstart);
            fprintf('epsilon=%.3f,t=%.4f/%.3f, dt=%.2e, Nx=%d, Ny=%d, timeElapsed=%f\n',epsilon,t,T,dt,Nx,Ny,timeElapsed);
        end
        
        if 1 == option.saveflag
            ss = [dir_data '/phi1_t=' num2str(t) '.txt'];
            fid1 = fopen(ss, 'wt');
            fprintf(fid1, '%f\n', phi(:));
            fclose(fid1);
        end
    end
        
end



if 1 == option.savefinal
    name=['phi_' pde.name '_e',num2str(pde.epsilon),'M',num2str(pde.M),...
        'Nx=',num2str(Nx),'Ny=',num2str(Ny),'dt=',num2str(dt)];
    filename=[name '_SAV_BDF.mat'];
    save(filename,'epsilon','xx','yy','hx','hy','Nx','Ny','dt','T','phi','domain');
end

end


function r = fun_r_init(phi)
global hx hy C0
E1 = fft2(F(phi));
if E1(1,1)*hx*hy + C0 <0
    fprintf("Root < 0, E + C=%f\n",E1(1,1)*hx*hy);
    return;
end
r  = sqrt(E1(1,1)*hx*hy+C0);
end

function r = fun_H(phi)
global hx hy C0
E1 = fft2(F(phi));
if E1(1,1)*hx*hy + C0 <0
    fprintf("Root < 0, E + C=%f\n",E1(1,1)*hx*hy);
    return;
end
r = fun(phi)./sqrt(E1(1,1)*hx*hy + C0);
end

function result = inv_A(phi)
global dt k2 M epsilon S
    L = -epsilon^2*k2 + S;
    phihat = fft2(phi);
    r      = phihat./(3/(2*M*dt) -k2 .* L);
    result = real(ifft2(r));
end

function r = fun(phi)
global S 
    r = phi.*(phi.^2 - 1 - S);
end

function r = F(phi)
global S
    r = (phi.^2 - 1).^2./4 - S/2*phi.^2;
end

function r = inner_int(f,g)
global hx hy
    E1 = fft2(f.*g);
    r  = E1(1,1)*hx*hy;
end

function r = diff_lap(phi)
global k2
    r=real(ifft2((k2.*fft2(phi))));
end

