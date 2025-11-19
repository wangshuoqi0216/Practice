function pde = ex01_CH_2D_data(para)

if nargin == 0
    epsilon = 1;
    M       = 1;
    alpha   = 1;
    beta = 1;
    C0 = 0;
    S  = 0;
else
    if ~isstruct(para)
        exit('we need a struct data');
    end
    if ~isfield(para,'epsilon') || isempty(para.epsilon)
        epsilon = 1;
    else
        epsilon = para.epsilon;
    end
    if ~isfield(para,'M') || isempty(para.M)
        M = 1;
    else
        M = para.M;
    end
    if ~isfield(para,'alpha') || isempty(para.alpha)
        alpha = 1;
    else
        alpha = para.alpha;
    end 
    if ~isfield(para,'beta') || isempty(para.beta)
        beta = 1;
    else
        beta = para.beta;
    end 
    if ~isfield(para,'C0') || isempty(para.C0)
        C0 = 1;
    else
        C0 = para.C0;
    end     
    if ~isfield(para,'S') || isempty(para.S)
        S = 1;
    else
        S = para.S;
    end     
end

pde = struct('epsilon',epsilon, ...
             'M', M, ...
             'alpha',alpha, ...
             'beta',beta, ...
             'C0',C0, ...
             'S', S, ...
             'initphi',@initphi, ...
             'name','ex01_CH_2D_data');

    function r = initphi(xx,yy)  
        theta = atan2(yy, xx);
        r = tanh((1.7 + 1.2*cos(6*theta) - sqrt(xx.^2 + yy.^2)) / (sqrt(2)*epsilon));
    end
end
%计算角度
%theta = atan2(yy, xx);
%初始条件：星形界面
%u = tanh((1.7 + 1.2*cos(6*theta) - sqrt(xx.^2 + yy.^2)) / (sqrt(2)*epsilon));