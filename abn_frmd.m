function  [Z, value, Lr] = ABN_FRMD(X,A, p,lambda,lambda0,rc,rho,para)


if nargin < 10
    display = true;
end
[~, n] = size(X);
XtX = p/2*(A'*A);
XtD = p/2*(A'*X);
maxiter = 200;
mu = rc*norm(X,2);
tol2 = 1e-3;
I = eye(n);
Z = ones(size(XtD));
Z_old = Z;
W1 = eye(n,n);
W2 = ones(n,1); 


[pairs,wcost,numpairs]=get_nn_graph(X,para.knn);
R = zeros(size(X,2),numpairs);
for i=1:numpairs
    R(pairs(1,i)+1,i) =  wcost(i);
    R(pairs(2,i)+1,i) = -wcost(i);
end
R = R/(para.knn-1);
L = 0.5*R*R';

Lr = L+para.elpson*eye(size(W1));


%
W = (lambda*Lr+p/2*lambda0*W1)*diag(W2);
for t = 1 : maxiter
   Z = lyap(XtX,W,-XtD);
   
   W1 = (Z'*Z+mu^2*I)^(p/2-1);


   E = X-A*Z;   
   E = dot(E,E);
   W2 =  (E+mu^2).^(1-p/2);   
   W = (lambda*Lr+p/2*lambda0*W1)*diag(W2);  
   
   % update mu
   mu = mu/rho; 
   
   value(t) = norm(Z_old-Z,'fro')/norm(Z,'fro');
   if value(t) <tol2       
       break;
   end 
   Z_old = Z;
end


