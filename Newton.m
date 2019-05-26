function y = Newton(x)
%as constantes usadas aqui ser�o para o tipo de solo arenoso
%par�metro do solo arenoso
b = 3.31;
%umidade do solo saturado
tetaS = 0.4;
%potencial do solo saturado
psiS = -0.093;
%largura m�dia da folha da maple tree
lf = 0.01;
%altura m�dia da maple tree
hc = 0.8;
%coeficiente de transfer�ncia 
Ct = 156.2;
%velocidade do vento m�dio (11,6 km/h - estimado)
Uhc = 41.76;
%index da area da folha
LAI = 3;
%coeficiente de arrasto (estimado)
cd = 1.28;

%uma vari�vel auxiliar para pegar o tamanho do vetor x
auxiliar = size(x);

%inicializando o vetor de sa�da y
y = repelem(0,2,1);

%uma possivel solu��o para o potencial da folha � este ser igual ao
%potencial da ra�z
y(1,1) = psiS*(x(auxiliar(1),1)/tetaS)^(-b);

%uma poss�vel solu��o para a solu��o da resist�ncia 

end