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
             'exactphi',@exactphi, ...
             'rhsphi',@rhsphi, ...
             'name','ex01_CH_2D_data');

    function r = exactphi(x,y,t)  
        r = sin(x).*sin(y).*exp(-t);
    end

    function r = initphi(x,y)  
        r = exactphi(x,y,0); 
    end

    function r = rhsphi(x,y,t)  
        r = -exp(-3.*t).*sin(x).*sin(y).*(exp(2.*t) + 2.*M.*exp(2.*t) + 6.*M.*sin(x).^2 + 6.*M.*sin(y).^2 - 18.*M.*sin(x).^2.*sin(y).^2 - 4.*M.*epsilon.^2.*exp(2.*t));
    end

end