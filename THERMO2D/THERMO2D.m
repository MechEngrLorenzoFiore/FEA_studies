% the program can solve 2D problems of thermal flow
clc,clear

% material thermal properties
kx = 10;
kz = 10;

% system geometric characteristics
hx = 1.25;
hz = 0.5;

% number of elements per axis
n_elem_x = 8;
n_elem_z = 2;
n_elem = n_elem_x * n_elem_z;





% ------------------------------ mesh ----------------------------------------
% matrix M contains the identificative number of the node in the mesh position
M(1,:) = [1 : 1 : n_elem_x+1];
for i = 2 : n_elem_z +1
    M(i,:) = M(i-1,:)+n_elem_x+1;
end
M = flip(M);

% number of nodes
n = max(max(M)); 





% --------------------- connections matrix -------------------------------
% matrix P contains, for every element, info about the elements connected
P = zeros(n_elem,5);
Maux = M(2:end,1:end-1);
Maux = flip(Maux);
Maux = Maux';
Maux = Maux(:);

for i=1:n_elem
    [r,c]=find(M==Maux(i));
    P(i,:) = [ i M(r,c) M(r,c+1) M(r-1,c+1) M(r-1,c)];
end






% ----------------------- loads and constrains -------------------------
% coloumn vector of every node in order
nn = [1:1:n]';

% row vector of the nodes at known temperature
n_Tv = [1 2 3 4 5 6 7 8 9 10 18 19 20 21 22 23 24 25 26 27];

% temperature is known in all the system boundary
Tv=[ 0 3.826 7.071 9.238 10 9.238 7.071 3.826 0 0 0 0 38.26 70.71 92.38 100 92.38 70.71 38.26 0 ];


% logic vectors identifing constrained and free positions
vin = (logical((sum((nn==n_Tv)'))'));
lib = ~vin;

% thermal flow vector inizialization
Ql = zeros( n-length(Tv),1);





% ------------------------- solution -------------------------------------
s_eta = [ -1 -1 +1 +1]';        
s_xi  = [ -1 +1 +1 -1]';        
s = [ s_eta s_eta s_xi s_xi];   

% Gauss point for numerical integration
val = 1/sqrt(3);
gp = [ +val -val +val -val ];   

fi=zeros(4,4,4);
for k=1:4
    aux = 0.5.*(ones(4,1) + gp(k) .* s(:,k) );
    fi(:,:,k) = aux * aux';
end

Dfi = zeros(4,4,2);
Dfi(:,:,1) = 0.25 * s_xi * s_xi';
Dfi(:,:,2) = 0.25 * s_eta * s_eta';

Jx = 0.5 * hx;
Jz = 0.5 * hz;




% ------------------------ stiffness matrix -----------------------------
K=zeros(n); 

ke(:,:) = 2*Jx*Jz*( kx/Jx^2 * ( fi(:,:,1) + fi(:,:,2) ) .* Dfi(:,:,1) +...
                    kz/Jz^2 * ( fi(:,:,3) + fi(:,:,4) ) .* Dfi(:,:,2) );

for i=1:n_elem  
    for j=1:4
        for f=1:4 
            K(P(i,1+j),P(i,1+f))=K(P(i,1+j),P(i,1+f))+ke(j,f);  
        end
    end   
end


Kll = K( lib , lib );
Klv = K( lib , vin );

Tl = Kll \ ( Ql - Klv*Tv' );

T = zeros(n,1);
T(lib) = Tl;
T(vin) = Tv;

