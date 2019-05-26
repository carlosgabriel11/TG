function y = Newton(x)
%as constantes usadas aqui serão para o tipo de solo arenoso
%parâmetro do solo arenoso
b = 3.31;
%umidade do solo saturado
tetaS = 0.4;
%potencial do solo saturado
psiS = -0.093;
%largura média da folha da maple tree
lf = 0.01;
%altura média da maple tree
hc = 0.8;
%coeficiente de transferência 
Ct = 156.2;
%velocidade do vento médio (11,6 km/h - estimado)
Uhc = 41.76;
%index da area da folha
LAI = 3;
%coeficiente de arrasto (estimado)
cd = 1.28;

%uma variável auxiliar para pegar o tamanho do vetor x
auxiliar = size(x);

%inicializando o vetor de saída y
y = repelem(0,2,1);

%uma possivel solução para o potencial da folha é este ser igual ao
%potencial da raíz
y(1,1) = psiS*(x(auxiliar(1),1)/tetaS)^(-b);

%uma possível solução para a solução da resistência 

end