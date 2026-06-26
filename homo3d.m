%%%% HOMOGENIZATION WITH SENSITIVITY ANALYSIS FOR TOPOLOGY OPTIMIZATION %%%%
function [DH,dDH] = homo3d(rho,lx,ly,lz,E0,Emin,nu,penal,edgeN)
%% MATERIAL INTERPOLATION
[ny, nx, nz] = size(rho); % Size of sub-cell mesh
E = Emin+rho.^penal*(E0-Emin); % SIMP
dE = penal*rho.^(penal-1)*(E0-Emin);
nelx = nx/edgeN; nely = ny/edgeN; nelz = nz/edgeN;
nele = nelx*nely*nelz;
%% GENERATE MESH & APPLY PERIODIC BOUNDARY CONDITIONS
[cornerNodes,edgeNodesX,edgeNodesY,edgeNodesZ] = getNodesPBC(nelx,nely,nelz);
enodeMat = zeros(nele,20);
for k=1:nelz
    for i=1:nelx
        for j=1:nely
             id = j+(i-1)*nely+(k-1)*nely*nelx; %element ID
             enodeMat(id,1) = cornerNodes(j+1,i,k);
             enodeMat(id,2) = cornerNodes(j+1,i+1,k);
             enodeMat(id,3) = cornerNodes(j,i+1,k);
             enodeMat(id,4) = cornerNodes(j,i,k);
             enodeMat(id,5) = cornerNodes(j+1,i,k+1);
             enodeMat(id,6) = cornerNodes(j+1,i+1,k+1);
             enodeMat(id,7) = cornerNodes(j,i+1,k+1);
             enodeMat(id,8) = cornerNodes(j,i,k+1);
             enodeMat(id,9) = edgeNodesX(j+1,i,k);
             enodeMat(id,10) = edgeNodesY(j,i+1,k);
             enodeMat(id,11) = edgeNodesX(j,i,k);
             enodeMat(id,12) = edgeNodesY(j,i,k);
             enodeMat(id,13) = edgeNodesZ(j+1,i,k);
             enodeMat(id,14) = edgeNodesZ(j+1,i+1,k);
             enodeMat(id,15) = edgeNodesZ(j,i+1,k);
             enodeMat(id,16) = edgeNodesZ(j,i,k);
             enodeMat(id,17) = edgeNodesX(j+1,i,k+1);
             enodeMat(id,18) = edgeNodesY(j,i+1,k+1);
             enodeMat(id,19) = edgeNodesX(j,i,k+1);
             enodeMat(id,20) = edgeNodesY(j,i,k+1);
        end
    end
end
%% DEGREES OF FREEDOM
edofMat = zeros(nelx*nely*nelz,60);
edofMat(:,1:3:end) = 3*enodeMat-2;
edofMat(:,2:3:end) = 3*enodeMat-1;
edofMat(:,3:3:end) = 3*enodeMat;
nnP = 4*nelx*nely*nelz;
ndof = 3*nnP; % Number of dofs
%% MAPPING SUB-CELLS TO PARENT ELEMENT INDICES (elabels) AND LOCAL SUB-CELL POSITIONS (plabels)
elabelsXY = kron(reshape(1:nelx*nely,nely,nelx),ones(edgeN,edgeN));
plabelsXY = repmat(reshape(1:edgeN*edgeN,edgeN,edgeN),nely,nelx);
elabels = zeros(ny,nx,nz);
plabels = zeros(ny,nx,nz);
for i=1:nelz
    for k=1:edgeN
        elabels(:,:,edgeN*(i-1)+k) = elabelsXY+(i-1)*nely*nelx;
        plabels(:,:,edgeN*(i-1)+k) = plabelsXY+(k-1)*edgeN*edgeN;
    end
end
%% ELEMENT MATRICES
dx = lx/nelx; dy = ly/nely; dz = lz/nelz;
[ke,keSub,fe,feSub] = stiffnessMatrix(dx,dy,dz,1,nu,edgeN);
%% GLOBAL MATRICES
nsub = nx*ny*nz;
iK = kron(edofMat,ones(60,1))';
jK = kron(edofMat,ones(1,60))';
iF = repmat(edofMat',6,1);
jF = [ones(60,nele); 2*ones(60,nele); 3*ones(60,nele); 4*ones(60,nele); 5*ones(60,nele); 6*ones(60,nele)];
kes = zeros(length(ke(:)),nele);
fes = zeros(length(fe(:)),nele);
for s = 1:nsub
    kes(:,elabels(s)) = kes(:,elabels(s))+keSub{plabels(s)}(:)*E(s);
    fes(:,elabels(s)) = fes(:,elabels(s))+feSub{plabels(s)}(:)*E(s);
end
sK = kes(:);
K = sparse(iK(:), jK(:), sK(:), ndof, ndof);K=(K+K')/2;
sF = fes(:);
F = sparse(iF(:), jF(:), sF(:), ndof, 6);
%% DISPLACEMENT FIELDS
chi(4:ndof,:) = K(4:ndof,4:ndof) \ F(4:ndof,:);
chi0 = zeros(nele, 60, 6);
chi0_e = zeros(60, 6);
chi0_e([4 7:11 13:60],:) = ke([4 7:11 13:60],[4 7:11 13:60])\fe([4 7:11 13:60],:);
chi0(:,:,1) = kron(chi0_e(:,1)', ones(nele,1));
chi0(:,:,2) = kron(chi0_e(:,2)', ones(nele,1));
chi0(:,:,3) = kron(chi0_e(:,3)', ones(nele,1));
chi0(:,:,4) = kron(chi0_e(:,4)', ones(nele,1));
chi0(:,:,5) = kron(chi0_e(:,5)', ones(nele,1));
chi0(:,:,6) = kron(chi0_e(:,6)', ones(nele,1));
%% HOMOGENIZATION
DH = zeros(6);
dDH = cell(ny,nx,nz);
dDH(:) = {zeros(6)};
cellVolume = lx*ly*lz;
for row =1:6
    lambda_row = (chi0(:,:,row)-chi(edofMat+(row-1)*ndof));
    for col = 1:6
        lambda_col = (chi0(:,:,col)-chi(edofMat+(col-1)*ndof));
        for s = 1:nsub
            val = (lambda_row(elabels(s),:)*keSub{plabels(s)})*lambda_col(elabels(s),:)';
            DH(row,col) = DH(row,col) + 1/cellVolume*E(s)*val;
            dDH{s}(row,col) = 1/cellVolume*dE(s)*val;
        end
    end
end
end

% THE FUNCTION IS USED TO OBTAIN THE NODE INFORMATION CONSIDERING PERIODIC
% BOUNDARY CONDITIONS
function [cornerNodes,edgeNodesX,edgeNodesY,edgeNodesZ] = getNodesPBC(nelx,nely,nelz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cornerNodes : Corner nodes
% edgeNodesX  : Nodes on X-edges
% edgeNodesY  : Nodes on Y-edges
% edgeNodesZ  : Nodes on Z-edges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
surfaceNodeNumUni = 3*nelx*nely;
cornerNodesFirstL = (0:2:2*(nely-1))'+(1:2*nely+1-1+nely+1-1:surfaceNodeNumUni);
cornerNodesFirstL(end+1,:) = cornerNodesFirstL(1,:);
cornerNodesFirstL(:,end+1) = cornerNodesFirstL(:,1);
edgeNodesXFirstL = (0:1:nely-1)'+(2*nely+1:3*nely:surfaceNodeNumUni);
edgeNodesXFirstL(end+1,:) = edgeNodesXFirstL(1,:);
edgeNodesYFirstL = (0:2:2*nely-2)'+(2:3*nely:surfaceNodeNumUni);
edgeNodesYFirstL(:,end+1) = edgeNodesYFirstL(:,1);
edgeNodesZSecondL = reshape(1:nely*nelx,nely,nelx) + surfaceNodeNumUni;
edgeNodesZSecondL(end+1,:) = edgeNodesZSecondL(1,:);
edgeNodesZSecondL(:,end+1) = edgeNodesZSecondL(:,1);
cornerNodes = zeros(nely+1,nelx+1,nelz+1);
edgeNodesX = zeros(nely+1,nelx,nelz+1);
edgeNodesY = zeros(nely,nelx+1,nelz+1);
edgeNodesZ = zeros(nely+1,nelx+1,nelz);
for i=1:nelz
    cornerNodes(:,:,i) = cornerNodesFirstL + (i-1)*(surfaceNodeNumUni+nelx*nely);%nnodes_z_layer_u=nelx*nely
    edgeNodesX(:,:,i) = edgeNodesXFirstL + (i-1)*(surfaceNodeNumUni+nelx*nely);
    edgeNodesY(:,:,i) = edgeNodesYFirstL + (i-1)*(surfaceNodeNumUni+nelx*nely);
end
cornerNodes(:,:,nelz+1) = cornerNodes(:,:,1);
edgeNodesX(:,:,nelz+1) = edgeNodesX(:,:,1);
edgeNodesY(:,:,nelz+1) = edgeNodesY(:,:,1);
for i=1:nelz
    edgeNodesZ(:,:,i) = edgeNodesZSecondL + (i-1)*(surfaceNodeNumUni+nelx*nely);%nnodes_z_layer_u=nelx*nely
end
end

% THE FUNCTION IS USED TO CALCULATE THE STIFFNESS MATRIX AND LOAD MATRIX
function [ke,keSub,fe,feSub] = stiffnessMatrix(dx,dy,dz,E0,nu,edgeN)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ke       : Stiffness matrix
% fe       : Load matrix
% keSub    : The stiffness matrix contribution of each sub-cell
% feSub    : The load matrix contribution of each sub-cell
% dx,dy,dz : The size of the element
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xyz = [0,0,0;dx,0,0;dx,dy,0;0,dy,0;0,0,dz;dx,0,dz;dx,dy,dz;0,dy,dz;
    dx/2,0,0;dx,dy/2,0;dx/2,dy,0;0,dy/2,0;0,0,dz/2;dx,0,dz/2;dx,dy,dz/2;0,dy,dz/2;
    dx/2,0,dz;dx,dy/2,dz;dx/2,dy,dz;0,dy/2,dz];
gauss = [-sqrt(3/5), 0, sqrt(3/5)];
weights = [5/9, 8/9, 5/9];
ke = zeros(60, 60);
keSub=cell(edgeN,edgeN,edgeN);
fe = zeros(60, 6);
feSub=cell(edgeN,edgeN,edgeN);
length_gauss = length(gauss);
D_left_top = [1 nu/(1-nu) nu/(1-nu);nu/(1-nu) 1 nu/(1-nu);nu/(1-nu) nu/(1-nu) 1];
D_right_bottom = diag(ones(3,1)*(1-2*nu)/2/(1-nu));
D = E0*(1-nu)/(1+nu)/(1-2*nu)*[D_left_top,zeros(3,3);zeros(3,3),D_right_bottom];
dx_sub = dx/(edgeN);
dy_sub = dy/(edgeN);
dz_sub = dz/(edgeN);
xis = [-1 1 1 -1 -1 1 1 -1 0 1 0 -1 -1 1 1 -1 0 1 0 -1];
etas = [-1 -1 1 1 -1 -1 1 1 -1 0 1 0 -1 -1 1 1 -1 0 1 0];
mus = [-1 -1 -1 -1 1 1 1 1 -1 -1 -1 -1 0 0 0 0 1 1 1 1];
for k = 1:edgeN
    for l = 1:edgeN
        for m = 1:edgeN
            x0 = (l-1) * dx_sub; %coordinate of the node at lower left corner
            y0 = (edgeN-k) * dy_sub;
            z0 = (m-1) * dz_sub;
            keSub{k,l,m}=zeros(60,60);
            feSub{k,l,m}=zeros(60,6);
            for ii = 1:length_gauss
                for jj = 1:length_gauss
                    for kk = 1:length_gauss
                        w = weights(ii)*weights(jj)*weights(kk);
                        xi_sub = gauss(ii);
                        eta_sub = gauss(jj);
                        mu_sub = gauss(kk);
                        x = (1-xi_sub)/2*x0+(1+xi_sub)/2*(x0+dx_sub);
                        y = (1-eta_sub)/2*y0+(1+eta_sub)/2*(y0+dy_sub);
                        z = (1-mu_sub)/2*z0+(1+mu_sub)/2*(z0+dz_sub);
                        xi = 2*x/dx-1;
                        eta = 2*y/dy-1;
                        mu = 2*z/dz-1;
                        dN_dxi = zeros(1,20);
                        dN_deta = zeros(1,20);
                        dN_dmu = zeros(1,20);
                        %% corner
                        for i=1:8
                            dN_dxi(i) = 0.125*(1+etas(i)*eta)*(1+mus(i)*mu)*(2*xi+xis(i)*etas(i)*eta+xis(i)*mus(i)*mu-xis(i));
                            dN_deta(i) = 0.125*(1+xis(i)*xi)*(1+mus(i)*mu)*(2*eta+etas(i)*xis(i)*xi+etas(i)*mus(i)*mu-etas(i));
                            dN_dmu(i) = 0.125*(1+xis(i)*xi)*(1+etas(i)*eta)*(2*mu+mus(i)*xis(i)*xi+mus(i)*etas(i)*eta-mus(i));
                        end
                        %% node 9,11,17,19
                        for i=[9,11,17,19]
                            dN_dxi(i) = -0.5*xi*(1+etas(i)*eta)*(1+mus(i)*mu);
                            dN_deta(i) = 0.25*(1-xi^2)*etas(i)*(1+mus(i)*mu);
                            dN_dmu(i) = 0.25*(1-xi^2)*(1+etas(i)*eta)*mus(i);
                        end
                        %% node 10,12,18,20
                        for i=[10,12,18,20]
                            dN_dxi(i) = 0.25*(1-eta^2)*xis(i)*(1+mus(i)*mu);
                            dN_deta(i) = -0.5*eta*(1+xis(i)*xi)*(1+mus(i)*mu);
                            dN_dmu(i) = 0.25*(1-eta^2)*(1+xis(i)*xi)*mus(i);
                        end
                        %% node 13,14,15,16
                        for i=13:16
                            dN_dxi(i) = 0.25*(1-mu^2)*xis(i)*(1+etas(i)*eta);
                            dN_deta(i) = 0.25*(1-mu^2)*etas(i)*(1+xis(i)*xi);
                            dN_dmu(i) = -0.5*mu*(1+xis(i)*xi)*(1+etas(i)*eta);
                        end
                        J = [dN_dxi;dN_deta;dN_dmu]*xyz;
                        detJ_sub = dx_sub*dy_sub*dz_sub/8;
                        dNdxyz = J \ [dN_dxi; dN_deta; dN_dmu];
                        dNdx = dNdxyz(1,:);
                        dNdy = dNdxyz(2,:);
                        dNdz = dNdxyz(3,:);
                        B = zeros(6,60);
                        for bi = 1:20
                            B(1,3*bi-2) = dNdx(bi);
                            B(2,3*bi-1)   = dNdy(bi);
                            B(3,3*bi)   = dNdz(bi);
                            B(4,3*bi-2) = dNdy(bi);
                            B(4,3*bi-1)   = dNdx(bi);
                            B(5,3*bi-1) = dNdz(bi);
                            B(5,3*bi)   = dNdy(bi);
                            B(6,3*bi-2) = dNdz(bi);
                            B(6,3*bi)   = dNdx(bi);
                        end
                        keSub{k,l,m} = keSub{k,l,m} + w * detJ_sub * B' * D * B;
                        feSub{k,l,m} = feSub{k,l,m} + w * detJ_sub * B' * D;
                    end
                end
            end
            ke=ke+keSub{k,l,m};
            fe=fe+feSub{k,l,m};
        end
    end
end
end

% ======================================================================= %
% A compact and efficient MATLAB code for homogenization method by using
% 20-node hexahedral elements with multiple sub-cells.
% Developed by: Chao Li et al.
% Email: chao.li-6@student.uts.edu.au
%
% Main reference:
% Dong, G., Tang, Y., Zhao, Y.F.: A 149 line homogenization code for 
% three-dimensional cellular materials written in matlab.
% Journal of Engineering Materials and Technology 141(1), 011005 (2019) 
%
% **************************   Disclaimer   ***************************** %
% The authors reserve all rights for the programs. The programs may be
% distributed and used for academic and educational purposes. The authors
% do not guarantee that the code is free from errors, and they shall not be
% liable in any event caused by the use of the program.
% ======================================================================= %
