close all;
clear; clc;

% add path
addpath('../','-begin');

% Space: Domain and N
N  = 64;

scheme_array = {'SAV_1st','SAV_CN','SAV_BDF'};

for ss = 1:length(scheme_array)
    scheme = scheme_array{ss};

    if 1 == strcmp(scheme,'SAV_1st')
        solver_fun = @solver_AllenCahn_2D_SAV_1st;
    elseif 1 == strcmp(scheme,'SAV_CN')
        solver_fun = @solver_AllenCahn_2D_SAV_CN;
    elseif 1 == strcmp(scheme,'SAV_BDF')
        solver_fun = @solver_AllenCahn_2D_SAV_BDF;
    end

    % PDE_array = {'data1','data2','data3'};
    PDE_array = {'data1'};


    for pp = 1:length(PDE_array)
        PDE = PDE_array{pp};

        % Parameters
        para.epsilon = 0.1;
        para.M = 1;
        para.S = 3;
        para.C0 = 40;
        para.name = PDE;

        if 1 == strcmp(PDE,'data1')
            domain.left   = 0;
            domain.right  = 2.*pi;
            domain.bottom = 0;
            domain.top    = 4.*pi;

            Nx = N; Ny = 2*N;

            T = 0.36;
            dt_array = 1e-4*[ 16 8 4 2 1 ]';
            dt_ref = dt_array(end)/10;
            pde = ex02_AC_2D_data(para);
        end

        t0 = 0;
        tsave = 0.2 * T;

        maxIt = length(dt_array);

        %% option
        option.plotflag   = 0;
        option.printflag  = 0;
        option.vtkflag    = 0;
        option.saveflag   = 0;
        option.savefinal  = 1;
        option.energyflag = 0;
        option.tol = 1e-14;
        option.tolit = 1e-11;
        option.maxit = 2000;

        %% Run:
        delete *.mat
        tic;
        if ~isfield(pde,'exactphi') || ~isfield(pde,'rhsphi')
            time = struct('T',T,'t0',t0,'dt',dt_ref,'tsave',tsave);
            solver_fun(pde,domain,Nx,Ny,time,option);
        end
        for k = 1:maxIt
            dt = dt_array(k);
            time = struct('T',T,'t0',t0,'dt',dt,'tsave',tsave);
            solver_fun(pde,domain,Nx,Ny,time,option);
        end
        toc;

        %% Compute order of convergence
        error=zeros(maxIt,1);   order  =zeros(maxIt,1);
        if ~isfield(pde,'exactphi') || ~isfield(pde,'rhsphi')
            name=['phi_' pde.name '_e',num2str(pde.epsilon),'M',num2str(pde.M),...
                'Nx=',num2str(Nx),'Ny=',num2str(Ny),'dt=',num2str(dt_ref)];
            filename=[name '_' scheme '.mat'];
            load(filename,'phi');
            phi_exact = phi;

            clear phi;
            fprintf('Accuracy test with reference solution.\n');
        else
            Lx = domain.right - domain.left;
            Ly = domain.top   - domain.bottom;
            hx = Lx/Nx;
            hy = Ly/Ny;
            x_c  = domain.left   + hx*(0:Nx-1);
            y_c  = domain.bottom + hy*(0:Ny-1);
            [xx, yy] = ndgrid(x_c,y_c);
            phi_exact= pde.exactphi(xx,yy,T);
            fprintf('Accuracy test with exact solution.\n');
        end
        for k = 1:maxIt
            dt = dt_array(k);
            name=['phi_' pde.name '_e',num2str(pde.epsilon),'M',num2str(pde.M),...
                'Nx=',num2str(Nx),'Ny=',num2str(Ny),'dt=',num2str(dt)];
            filename=[name '_' scheme '.mat'];
            load(filename,'phi','hx','hy');
            err = fft2((phi_exact - phi).^2);
            error(k,1) = sqrt(err(1,1)*hx*hy);   % L2

            clear phi;
        end
        order(2:maxIt,:) = log(error(1:maxIt-1,:)./error(2:maxIt,:))./log(dt_array(1:maxIt-1,:)./dt_array(2:maxIt,:));

        %% Display error and order
        fprintf([func2str(solver_fun), '\n']);
        fprintf('%s  (epsilon=%g, S=%d, M=%.2e, Nx=%d, Ny=%d, T=%.2f)\n', pde.name, pde.epsilon, pde.S, pde.M, Nx, Ny, T);
        fprintf('   dt      &  Err_phi   & Order \n');
        for k = 1:maxIt
            fprintf('%.4e   %.4e   %.2f\n',dt_array(k),error(k,1),order(k,1));
        end
        fprintf('\n')

        % %% Plot
        % figure(2)
        % hh=loglog(dt_array,10*dt_array.^2,'k:','LineWidth',2);
        % xlabel('Time step $\delta t$','Interpreter','latex');
        % ylabel('$L^2$ error','Interpreter','latex');
        % grid on;
        % hold on;
        % hh=loglog(dt_array,error(:,1),'g>-');

        % %% Save error and order
        % name=['phi_e',num2str(para.epsilon),'M',num2str(para.M),...
        %       'Nx=',num2str(N),'Ny=',num2str(N)];
        % fileID = fopen([name,'.txt'],'w');
        % % fprintf(fileID,'%6s\n','%% Results');
        % % fprintf(fileID,'%6s\n','% dt	   &   Error_L2	   &  Order');
        % % A = [dt_array error];
        % % fprintf(fileID,'%.12f   %.4e   \n',A');
        % fprintf(fileID,'%.12e     %.4e      %.2f \n',[dt_array,error,order]');
        % fclose(fileID);
    end
end

%% results：
% Elapsed time is 50.136834 seconds.
% Accuracy test with reference solution.
% solver_AllenCahn_2D_SAV_1st
% ex02_AC_2D_data  (epsilon=0.1, S=3, M=1.00e+00, Nx=64, Ny=128, T=0.36)
%    dt      &  Err_phi   & Order 
% 1.6000e-03   1.0450e-02   0.00
% 8.0000e-04   5.2424e-03   1.00
% 4.0000e-04   2.6006e-03   1.01
% 2.0000e-04   1.2701e-03   1.03
% 1.0000e-04   6.0235e-04   1.08
% 
% Elapsed time is 76.955389 seconds.
% Accuracy test with reference solution.
% solver_AllenCahn_2D_SAV_CN
% ex02_AC_2D_data  (epsilon=0.1, S=3, M=1.00e+00, Nx=64, Ny=128, T=0.36)
%    dt      &  Err_phi   & Order 
% 1.6000e-03   3.3628e-05   0.00
% 8.0000e-04   8.4244e-06   2.00
% 4.0000e-04   2.1075e-06   2.00
% 2.0000e-04   5.2619e-07   2.00
% 1.0000e-04   1.3061e-07   2.01
% 
% Elapsed time is 62.159072 seconds.
% Accuracy test with reference solution.
% solver_AllenCahn_2D_SAV_BDF
% ex02_AC_2D_data  (epsilon=0.1, S=3, M=1.00e+00, Nx=64, Ny=128, T=0.36)
%    dt      &  Err_phi   & Order 
% 1.6000e-03   5.0232e-05   0.00
% 8.0000e-04   1.2560e-05   2.00
% 4.0000e-04   3.1389e-06   2.00
% 2.0000e-04   7.8330e-07   2.00
% 1.0000e-04   1.9437e-07   2.01


