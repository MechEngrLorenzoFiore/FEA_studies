% the program takes as input in the matrix P:
    % number of elements
    % angle of each beam (in degrees)
    % length
    % identificative number of the degree of freedom on the left
    % identificative number of the degree of freedom on the rigth
% and returns as output the stiffness matrix and the solution of
% the system in term of forces at the constained nodes

clc,clear % clean MATLAB workspace














% ---------------------------------------------------------------------
% -------------------- PRE PROCESSOR ----------------------------------
% input regarding geometry of the system
P=[1   90    L/(2*sqrt(2))    1   2    3   4;
   2   0     L/(2*sqrt(2))    5   6    1   2;
   3   45    L/2              5   6    3   4;
   4   -45   L/2              7   8    3   4;
   5   90    L*sqrt(2)/2      5   6    7   8;
   6   0     L*sqrt(2)/2      9   10   5   6;
   7   -45   L                11  12   5   6;
   8   45    L                9   10   7   8;
   9   0     L/sqrt(2)        11  12   7   8;
   10  -45   L                13  14   7   8];

% mechanical properties of the system
L = 2000;               % [mm] 
E = 217000;             % [N/mm^2] Young modulus
A = pi/4*(34^2-31^2);   % [mm^2]   beam section
% -------------------- END OF PRE PROCESSOR ---------------------------
% ---------------------------------------------------------------------














% ---------------------------------------------------------------------
% -------------------- SOLUTION ---------------------------------------

[s,~]=size(P);
n=max(max(P(:,[4:7])));   

% matrix initialization
K=zeros(n);     % stiffness matrix
k=zeros(4,4,s); % auxiliar matrix 

% main loop to build the stiffness matrix from the known informations
for i = 1:s  
    alp = P(i,2)*2*pi/360;
    Lt = P(i,3);   
    k(:,:,i)=E*A/Lt*[cos(alp)^2  cos(alp)*sin(alp)  -cos(alp)^2  -cos(alp)*sin(alp);
             cos(alp)*sin(alp)  sin(alp)^2 -cos(alp)*sin(alp) -sin(alp)^2;
             -cos(alp)^2  -cos(alp)*sin(alp)  cos(alp)^2  cos(alp)*sin(alp);
             -cos(alp)*sin(alp)  -sin(alp)^2  cos(alp)*sin(alp) sin(alp)^2];
    for j=1:4
        for f=1:4  
            K(P(i,3+j),P(i,3+f))=K(P(i,3+j),P(i,3+f))+k(j,f,i);  
        end
    end   
end

% extract the characteristic matrices for free and constrained nodes from K
KLL=K([1:8],[1:8]);
KLV=K([1:8],[9:14]);
KVL=K([9:14],[1:8]);
KVV=K([9:14],[9:14]);

% known forces and displacements
FL=[26000 -72000 0 0 0 0 0 0]';
dV=zeros(6,1);

% solution for displacements at free nodes
dL = KLL^(-1)*(FL-KLV*dV);

% solution for reactions at constrained nodes
FV = KVL*dL+KVV*dV;

% -------------------- END OF SOLUTION --------------------------------
% ---------------------------------------------------------------------














% ---------------------------------------------------------------------
% -------------------- POST PROCESSOR ---------------------------------
d = [dL; dV];    % vector of all the displacements
q = zeros(4,1);  % auxiliary vector of the displacements of the nodes of a single beam
Q = zeros(4,s);  % auxiliary vector of the displacements of all the nodes 

for i=1:s       
    for j=1:4  
        q(j)=d(P(i,3+j));
    end
    Q(:,i)=k(:,:,i)*q;
end

% the sign of the force component on the left node of the beam 
% is the sign of the stress on it (tension positive, compression negative)
R = 10^(-3).*Q([3:4],:)';

% force on the beams
R(:,3) = (R(:,1).^2 + R(:,2).^2).^(1/2);

% stress on the beams
R(:,4)=10^(3)*R(:,3)./A;

% --------------------END OF POST PROCESSOR ---------------------------
% ---------------------------------------------------------------------


















%benchmark:
% L=1;
% E=1;
% A=1;
% P=[1  0  L  1  2   3  4; 
%    2  0  L  3  4   5  6];
% P=[1  0     L          5  6   1  2; 
%    2  45    L*sqrt(2)  3  4   1  2;
%    3  -45   L*sqrt(2)  7  8   1  2;
%    ];


%     brutalmente:
%     K([P(i,4):P(i,5)],[P(i,4):P(i,5)])=K([P(i,4):P(i,5)],[P(i,4):P(i,5)])+k([1:2],[1:2]);
%     K([P(i,4):P(i,5)],[P(i,6):P(i,7)])=K([P(i,4):P(i,5)],[P(i,6):P(i,7)])+k([1:2],[3:4]);
%     K([P(i,6):P(i,7)],[P(i,4):P(i,5)])=K([P(i,6):P(i,7)],[P(i,4):P(i,5)])+k([3:4],[1:2]);
%     K([P(i,6):P(i,7)],[P(i,6):P(i,7)])=K([P(i,6):P(i,7)],[P(i,6):P(i,7)])+k([3:4],[3:4]);
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
