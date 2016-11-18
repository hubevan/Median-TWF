%% Implementation of the median truncated Wirtinger Flow (median-TWF) algorithm which is adapted from 
%  TWF by Y. Chen and E. J. Candès.

function [Relerrs] = medianTWF(y, x, Params, A, At)    
%% Initialization
    Relerrs=zeros(Params.T+1,1);
    npower_iter = Params.npower_iter;           % Number of power iterations 
    z0 = randn(Params.n1,Params.n2); z0 = z0/norm(z0,'fro');    % Initial guess 
    normest = sqrt(median(y(:))/0.455);    % Estimate norm to scale eigenvector  
    
    for tt = 1:npower_iter,                     % Truncated power iterations
        ytr = y.* (abs(y) <= Params.alpha_y^2 * normest^2 );
        z0 = At( ytr.* A(z0) ); z0 = z0/norm(z0,'fro');
    end
    
    z = normest * z0;                   % Apply scaling 
    Relerrs(1) = norm(x - exp(-1i*angle(trace(x'*z))) * z, 'fro')/norm(x,'fro'); % Initial rel. error
    
    %% Loop
    grad_type = Params.grad_type;
    if strcmp(grad_type, 'TWF_Poiss') == 1
        mu = @(t) Params.mu; % Schedule for step size 
    elseif strcmp(grad_type, 'WF_Poiss') == 1
        tau0 = 330;                         % Time constant for step size
        mu = @(t) min(1-exp(-t/tau0), 0.2); % Schedule for step size  
    end
    m=Params.m;
    for t = 1: Params.T,
        yz=A(z);
        absyz=abs(yz);
        Kt=median(abs(absyz.^2-y(:)))/0.4;

        if strcmp(Params.grad_type,'TWF_Poiss') == 1   % truncated gradient / Wirtinger flow
            % truncation rules
            Eub =  abs(yz) / norm(z)   <= Params.alpha_ub;
            Elb =  abs(yz) / norm(z)   >= Params.alpha_lb;
            Eh  =  abs(y - abs(yz).^2) <= Params.alpha_h * Kt / norm(z) * absyz;
            %Eh  =  abs(y - abs(yz).^2) <= 3 * Kt / norm(z) * abs(yz);

            grad  = 1/m* At( 2* ( absyz.^2-y ) ./ (absyz.^2) .*yz ...
                              .* Eub .* Elb .* Eh );    % truncated Poisson gradient

        elseif strcmp(Params.grad_type,'WF_Poiss') == 1    % untruncated gradient / Wirtinger flow
            grad  = 1/m* At( 2* ( absyz.^2-y ) ./ (absyz.^2) .*yz ); % Poisson gradient
        end
        
        z = z - mu(t) * grad;             % Gradient update 
        
        currenterrs=norm(x - exp(-1i*angle(trace(x'*z))) * z, 'fro')/norm(x,'fro');
        
        if currenterrs<1e-10
            break;
        end
        Relerrs(t+1) = currenterrs;  
    end
    for t=t:Params.T,
        Relerrs(t+1)= currenterrs;  
    end