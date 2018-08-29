g=[500 300];
f=-1*g;
A=[20 10;
   3 2];
b=[1000;   180];
[x,z]=linprog(f,A,b);

fprintf('En una granja agricola se desea criar conejos y pollos como complemento\n en su economia de forma que no se superen en conjunto las 180 horas\n mensuales destinadas a esta actividad. Su almacén sólo puede albergar un\n máximo de 1000 kilogramos de pienso. Si se supone que un conejo necesita\n 20 kilogramos de pienso al mes y un pollo 10 kilogramos al mes que las horas\n mensuales de cuidados requeridos por un conejo son 3 y por un pollo son\n 2 y que los beneficios que reportaría su venta ascienden a 500 y 300 pesetas\n por cabeza respectivamente hallar el número de animales que deben\n criarse para que el beneficio sea máximo\n');
fprintf('\n El vector a maximizar es:\n');
disp(g);
fprintf('\n La matriz de restricciones es: \n');
disp(A);
fprintf('\n La matriz de resultados es: \n');
disp(b);
fprintf('\n Los valores optimos son los siguientes\n');
disp(x);
fprintf('\n Para conseguir: ');
disp(-z);
