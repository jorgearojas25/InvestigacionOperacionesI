f=[60 80];
A=[ 1 0;
    0 1;
    1 1;
    -40 -50];
b=[ 8;
    10;
    9;
    -400];
[x,z]=linprog(f,A,b);

fprintf('los 400 alumnos de un colegio van a ir de excursi�n para ello se contrata el viaje a una empresa que dispone de 8autobuses del tipo A con 40 plazas y 10 del tipo B con 50 plazas pero \n solo de 9 conductores para ese d�a. Dada la diferente capacidad y calidad \n el alquiler de cada autob�s de los grandes (tipo B) cuesta 80 � y el de cada uno \n de los peque�os (tipo A) 60 � �Cu�ntos autobuses de \n cada clase convendr� alquilar para que el viaje resulte lo m�s econ�mico posible?');
fprintf('\n El vector a maximizar es:\n');
disp(f);
fprintf('\n La matriz de restricciones es: \n');
disp(A);
fprintf('\n La matriz de resultados es: \n');
disp(b);
fprintf('\n Los valores optimos son los siguientes\n');
disp(x);
fprintf('\n Para conseguir: ');
disp(z);
