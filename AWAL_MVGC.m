function [UU,A,Z,iter,obj,alpha,U] = AWAL_MVGC(X,Y,~,anchor,lambda1,lambda2)
% m      : the number of anchor. the size of Z is m*n.
% lambda1: the hyper-parameter of anchor weight regularization term.
% lambda2 : the hyper-parameter of structure alignment regularization term.
%anchor : m*d

%% initialize
maxIter = 50 ; % the number of iterations
numclass = length(unique(Y));
numview = length(X);
numsample = size(Y,1);
W = cell(numview, 1);%% Weight initialization
M_corr = cell(numview, 1); % Anchor correlation matrix
M_aux = cell(numview,1);
p = zeros(numclass,numsample);
p(:,1:numclass) = eye (numclass);
[Z, ~, ~]=svd(p','econ');
Z = Z';

for i = 1:numview
   [m,~] = size(anchor{i});
   W{i} = eye(m); 
   A{i} = anchor{i};
   A{i} = A{i}';
   X{i} = mapstd(X{i}',0,1); % turn into d*n
   M_aux{i} = (A{i}'*X{i}); 
   M_corr{i} = (A{i}' * A{i}); 
   for ii=1:numsample
       idx = 1:m;
       pp(idx,ii) = EProjSimplex_new(M_aux{i} (idx,ii)); 
   end
   Zi{i} = pp;
   clear pp;
end

alpha = ones(1,numview)/numview;
opt.disp = 0;

flag = 1;
iter = 0;
%% alternative optimization
while flag
    iter = iter + 1;
    
     %% optimize Ri
    for j=1:numview
        [Unew1,~,Vnew1] = svd(Z*Zi{j}','econ');
        R{j} = Unew1*Vnew1';
    end
    
     %% optimize W (weighted anchors)
    for iv = 1:numview
        % Quadratic programming to optimize w_diag
        [n,~] = size(anchor{iv});
        G = alpha(iv)*(Zi{iv} * Zi{iv}') .* M_corr{iv} + lambda1 * M_corr{iv}; 
        f = -2 * alpha(iv)* diag(Zi{iv}*X{iv}'*A{iv});
        options = optimoptions('quadprog', 'Display', 'off');
        w_diag = quadprog(G, f, [], [], ones(1, n), n, zeros(n, 1), [],[],options);
        % Update W{iv} as a diagonal matrix
        W{iv} = diag(w_diag);
    end

    %% optimize Z-i with weighted anchors
    for a = 1:numview
        % Weighted anchor matrix: A{a} * W{a}
        weighted_A = A{a} * W{a};
        M_aux{a} = (alpha(a)* weighted_A' * X{a} + lambda2 * R{a}' * Z) / (alpha(a) + lambda2); 
        for ii=1:numsample
            [n,~] = size(anchor{a});
            idx = 1:n;
            pp(idx,ii) = EProjSimplex_new(M_aux{a} (idx,ii));
        end
        Zi{a} = pp;
        clear pp;
    end

    %% optimize Z
    part1 = 0;
    for ia = 1:numview
        part1 = part1 + R{ia} *Zi{ia}; 
    end
    [Unew,~,Vnew] = svd(part1,'econ');
    Z = Unew*Vnew';

    %% optimize alpha
    P = zeros(numview,1);

    for iv = 1:numview
        % Weighted reconstruction error
        weighted_A = A{iv} * W{iv};
        P(iv) = norm( X{iv} - weighted_A * Zi{iv},'fro')^2;
    end
    Mfra = P.^-1;
    Q = 1/sum(Mfra);
    alpha = Q*Mfra;
    %%
    term1 = 0;
    term2 = 0;
    term3 = 0;
    for iv = 1:numview
        % Weighted reconstruction error
        weighted_A = A{iv} * W{iv};
        term1 = term1 + alpha(iv) * norm(X{iv} - weighted_A * Zi{iv},'fro')^2;
        term2 = term2 + lambda2 * norm(R{iv} *Zi{iv} - Z,'fro')^2; 
        % Weight regularization
        term3 = term3 + lambda1 * trace(W{iv}' * M_corr{iv} * W{iv});
    end
    obj(iter) = term1+ term2+term3;
    U{iter}=Z';
    if (iter>9) && (abs((obj(iter-1)-obj(iter))/(obj(iter-1)))<1e-6 || iter>maxIter || obj(iter) < 1e-10) 
        UU = Z'; 
        flag = 0;
    end
end
         
         
    
