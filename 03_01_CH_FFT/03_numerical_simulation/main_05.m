close all;
clear; clc;

% add path
addpath('../','-begin');

% Space: Domain and N
N  = 64;

solver_fun = @solver_CahnHilliard_2D_SAV_BDF;



    % PDE_array = {'data1','data2','data3'};
 PDE_array = {'data1'};


    for pp = 1:length(PDE_array)
        PDE = PDE_array{pp};

        % Parameters
        para.epsilon = 0.1;
        para.M = 1e-2;
        para.S = 3;
        para.C0 = 20;
        para.name = PDE;

        if 1 == strcmp(PDE,'data1')
            domain.left   = -pi;
            domain.right  = pi;
            domain.bottom = -pi;
            domain.top    = pi;

            Nx = 2*N; Ny = 2*N;


            T = 1;
            
            dt = 1e-3;
            pde = ex01_CH_2D_data(para);
        end

        t0 = 0;
        tsave = 0.2 * T;

    end

        %% option
        option.plotflag   = 0;
        option.printflag  = 0;
        option.vtkflag    = 0;
        option.saveflag   = 1;
        option.savefinal  = 1;
        option.energyflag = 0;
        option.tol = 1e-14;
        option.tolit = 1e-11;
        option.maxit = 2000;

       
       tic;
        
           
       time = struct('T',T,'t0',t0,'dt',dt,'tsave',tsave);
       solver_fun(pde,domain,Nx,Ny,time,option);
      
       toc;
