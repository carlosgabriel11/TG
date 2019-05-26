clear
clc

%parametro de cada solo
%aqui sera usado como exemplo o solo arenoso
b = 3.31;

%saturação da água no solo arenoso
tetaSat = 0.4;

%teta utilizado para linearização
tetaLin = 0.1;

%saturação do potencial de água no solo arenoso
psiSat = -0.093;

%condutividade de água no solo saturado de água e arenoso
kSat = 9.39*10^(-6);

%permeabilidade da raíz
kr = 10^(-8);

%raio da raíz
r_r = 10^(-4);

%número de discretizações verticais
nz = 10;

%volume de controle
R = 0.001: (0.01-0.001)/nz: 0.01 - (0.01-0.001)/nz;


%calculando algumas constantes
c1 = kSat*psiSat/((-1 - 3/b)*tetaLin^(b+3));

c3 = c1*(b+3)*((tetaLin)^(b + 2));

c5 = kSat*(2*b + 3)/tetaLin;

c6 = psiSat*(1+b);

c7 = -b*psiSat/tetaLin;

%número de discretizações na direção radial
nr = 10;

%discretização da raiz
r0 = (R)/nr;

%tamanho médio de profundidade da raíz
Z = 1;

%discretização vertical da raíz
z0 = Z/nz;

%criando a matriz Z
Z = zeros(nz*nr + 1, 1);

%criando a matriz A
A = zeros(nz*nr + 1);

%criando a matriz B
B = zeros(nr*nz + 1, 1);
B(2,1) = 1;

%criando a matriz C
C = repelem(0,1, nz*nr + 1);
C(1,nr*nz+1) = 1;

%elemento de passagem de água
nd = 5;

%preenchendo a matriz A na primeira camada
%primeiro elemento radial
A(1,1) = A(1,1)-2*c3/r0(1)^2-c3/z0^2+c5/(2*z0);
A(1,2) = A(1,2)+c3/(2*r0(1)^2) + c3/(r0(1)^2);
A(1, nz*nr + 1) = A(1, nz*nr + 1) - c3/(2*r0(1)^2) + c3/(r0(1)^2);
A(nz*nr + 1, 1) = A(nz*nr + 1, 1) + c3/(2*r_r^2) + c3/(r_r^2);
A(1,nr+1) = A(1,nr+1) +c3/z0^2-c5/(2*z0);
%Z(nz*nr + 1,1) = Z(nz*nr + 1,1) + 2*(-kr*z0 - kr*c6)/(pi*r_r);

%demais elementos radiais exceto o último
for counterRadius = 2:(nr-1)
    A(counterRadius, counterRadius - 1) = A(counterRadius, counterRadius - 1)-c3/(2*counterRadius*r0(1)^2) + c3/(r0(1)^2);
    A(counterRadius, counterRadius) = A(counterRadius, counterRadius)-2*c3/r0(1)^2;
    A(counterRadius, counterRadius + 1) = A(counterRadius, counterRadius + 1)+c3/(2*counterRadius*r0(1)^2) + c3/(r0(1)^2);
    
    %troca vertical
A(nr, nr) = A(nr, nr)-2*c3/r0(1)^2 + c3/(2*nr*r0(1)^2) + c3/(r0(1)^2);
    A(counterRadius, counterRadius+nr) = A(counterRadius, counterRadius+nr) +c3/z0^2-c5/(2*z0);
    A(counterRadius, counterRadius) = A(counterRadius, counterRadius)-c3/z0^2+c5/(2*z0);
    
end

%último elemento radial da primeira camada
A(nr,nr - 1) = A(nr,nr - 1)-c3/(2*nr*r0(1)^2) + c3/(r0(1)^2);
A(nr, nr) = A(nr, nr)-2*c3/r0(1)^2 + c3/(2*nr*r0(1)^2) + c3/(r0(1)^2) - c3/z0^2+c5/(2*z0);
A(nr, nr+nr) = A(nr, nr+nr) + c3/z0^2-c5/(2*z0);


%demais camadas verticais
for counterHeight = 2:(nz - 1)
    %primeira camada radial
    A((counterHeight - 1)*nr + 1, (counterHeight - 1)*nr + 1) = A((counterHeight - 1)*nr + 1, (counterHeight - 1)*nr + 1) - 2*c3/r0(counterHeight)^2 -c3*2/z0^2;
    A((counterHeight - 1)*nr + 1, (counterHeight - 1)*nr + 2) = A((counterHeight - 1)*nr + 1, (counterHeight - 1)*nr + 2) + c3/(2*r0(counterHeight)^2) + c3/(r0(counterHeight)^2);
    A((counterHeight - 1)*nr + 1, nr*nz + 1) = A((counterHeight - 1)*nr + 1, nr*nz + 1)  - c3/(2*r0(counterHeight)^2) + c3/(r0(counterHeight)^2);
    A(nz*nr + 1, (counterHeight - 1)*nr + 1) = A(nz*nr + 1, (counterHeight - 1)*nr + 1) + c3/(2*r_r^2) + c3/(r_r^2);
    %Z(nz*nr + 1,1) = Z(nz*nr + 1,1) + 2*(-kr*counterHeight*z0 - kr*c6)/(pi*r_r);
    A((counterHeight - 1)*nr + 1, (counterHeight - 1)*nr + 1 - nr) = A((counterHeight - 1)*nr + 1, (counterHeight - 1)*nr + 1 - nr) + c3/z0^2 + c5/(2*z0);
    A((counterHeight - 1)*nr + 1, (counterHeight - 1)*nr + 1 + nr) = A((counterHeight - 1)*nr + 1, (counterHeight - 1)*nr + 1 + nr) + c3/z0^2 - c5/(2*z0);
    
    %demais camadas radiais
    for counterRadius = 2:(nr - 1)
        A((counterHeight - 1)*nr + counterRadius, (counterHeight - 1)*nr + counterRadius - 1) = A((counterHeight - 1)*nr + counterRadius, (counterHeight - 1)*nr + counterRadius - 1)-c3/(2*counterRadius*r0(counterHeight)^2) + c3/(r0(counterHeight)^2);
        A((counterHeight - 1)*nr + counterRadius, (counterHeight - 1)*nr + counterRadius) = A((counterHeight - 1)*nr + counterRadius, (counterHeight - 1)*nr + counterRadius)-2*c3/r0(counterHeight)^2 - 2*c3/z0^2;
        A((counterHeight - 1)*nr + counterRadius, (counterHeight - 1)*nr + counterRadius + 1) = A((counterHeight - 1)*nr + counterRadius, (counterHeight - 1)*nr + counterRadius + 1)+c3/(2*counterRadius*r0(counterHeight)^2) + c3/(r0(counterHeight)^2);
        A((counterHeight - 1)*nr + counterRadius, (counterHeight - 1)*nr + counterRadius - nr) = A((counterHeight - 1)*nr + counterRadius, (counterHeight - 1)*nr + counterRadius - nr) + c3/z0^2 + c5/(2*z0);
        A((counterHeight - 1)*nr + counterRadius, (counterHeight - 1)*nr + counterRadius + nr) = A((counterHeight - 1)*nr + counterRadius, (counterHeight - 1)*nr + counterRadius + nr) + c3/z0^2 - c5/(2*z0);
        
    end
    
    %última camada radial
    A((counterHeight - 1)*nr + nr, (counterHeight - 1)*nr + nr) = A((counterHeight - 1)*nr + nr, (counterHeight - 1)*nr + nr)-2*c3/r0(counterHeight)^2 + c3/(2*nr*r0(counterHeight)^2) + c3/(r0(counterHeight)^2) - 2*c3/z0^2;
    A((counterHeight - 1)*nr + nr, (counterHeight - 1)*nr + nr - 1) = A((counterHeight - 1)*nr + nr, (counterHeight - 1)*nr + nr - 1)-c3/(2*nr*r0(counterHeight)^2) + c3/(r0(counterHeight)^2);
    A((counterHeight - 1)*nr + nr, (counterHeight - 1)*nr + nr - nr) = A((counterHeight - 1)*nr + nr, (counterHeight - 1)*nr + nr - nr) +c3/z0^2 + c5/(2*z0);
    A((counterHeight - 1)*nr + nr, (counterHeight - 1)*nr + nr + nr) = A((counterHeight - 1)*nr + nr, (counterHeight - 1)*nr + nr + nr) +c3/z0^2 - c5/(2*z0);

end

%última camada vertical
%primeira camada radial
A((nz - 1)*nr + 1, (nz - 1)*nr + 1) = A((nz - 1)*nr + 1, (nz - 1)*nr + 1) - 2*c3/r0(nz)^2 -2*c3/z0^2;% -c5/(2*z0) + c3/z0^2;
A((nz - 1)*nr + 1, (nz - 1)*nr + 2) = A((nz - 1)*nr + 1, (nz - 1)*nr + 2)+c3/(2*r0(nz)^2) + c3/(r0(nz)^2);
A((nz - 1)*nr + 1, nr*nz + 1) = A((nz - 1)*nr + 1, nr*nz + 1) - c3/(2*r0(nz)^2) + c3/(r0(nz)^2);
A(nz*nr + 1, (nz - 1)*nr + 1) = A(nz*nr + 1, (nz - 1)*nr + 1) + c3/(2*r_r^2) + c3/(r_r^2);
A((nz - 1)*nr + 1, (nz - 1)*nr + 1 - nr) = A((nz - 1)*nr + 1, (nz - 1)*nr + 1-nr)+c5/(2*z0) + c3/z0^2;
%Z(nz*nr + 1,1) = Z(nz*nr + 1,1) + 2*(-kr*nz*z0 - kr*c6)/(pi*r_r);

%demais elementos radiais exceto o último
for counterRadius = 2:(nr-1)
    A((nz - 1)*nr + counterRadius, (nz - 1)*nr + counterRadius - 1) = A((nz - 1)*nr + counterRadius, (nz - 1)*nr + counterRadius - 1)-c3/(2*counterRadius*r0(nz)^2) + c3/(r0(nz)^2);
    A((nz - 1)*nr + counterRadius, (nz - 1)*nr + counterRadius) = A((nz - 1)*nr + counterRadius, (nz - 1)*nr + counterRadius)-2*c3/r0(nz)^2 -2*c3/z0^2;
    A((nz - 1)*nr + counterRadius, (nz - 1)*nr + counterRadius + 1) = A((nz - 1)*nr + counterRadius, (nz - 1)*nr + counterRadius + 1)+c3/(2*counterRadius*r0(nz)^2) + c3/(r0(nz)^2);
    A((nz - 1)*nr + counterRadius, (nz - 1)*nr + counterRadius - nr) = A((nz - 1)*nr + counterRadius, (nz - 1)*nr + counterRadius - nr) +c5/(2*z0) + c3/z0^2;  
end

%última camada radial
A((nz - 1)*nr + nr, (nz - 1)*nr + nr) = A((nz - 1)*nr + nr, (nz - 1)*nr + nr)-2*c3/r0(nz)^2 + c3/(2*nr*r0(nz)^2) + c3/(r0(nz)^2) -2*c3/z0^2 ;
A((nz - 1)*nr + nr, (nz - 1)*nr + nr - 1) = A((nz - 1)*nr + nr, (nz - 1)*nr + nr - 1) -c3/(2*nr*r0(nz)^2) + c3/(r0(nz)^2);
A((nz - 1)*nr + nr, (nz - 1)*nr + nr - nr) = A((nz - 1)*nr + nr, (nz - 1)*nr + nr - nr) +c3/z0^2 +c5/(2*z0);

%sobre a raíz
A(nz*nr + 1, nz*nr + 1) = A(nz*nr + 1, nz*nr + 1) -2*c3/r_r^2 - c3/(2*r_r^2) + c3/(r_r^2);
A(nz*nr + 1, nz*nr + 1) = nz*A(nz*nr + 1, nz*nr + 1);