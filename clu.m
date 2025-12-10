function [idx] = clu(L,K)

L = (L + L')/2;
D = diag(1./sqrt(sum(L,2)));
L = D*L*D;
[U,S,V] = svd(L);

V = U(:,1:K);
V = D*V;

V = real(V);            
V = double(V);           
V(~isfinite(V)) = 0;      


idx = kmeans(V,K,'emptyaction','singleton','replicates',10,'display','off');
idx = idx';