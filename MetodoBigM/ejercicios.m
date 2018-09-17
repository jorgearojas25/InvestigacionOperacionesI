clear all;
clc;


%% Ejemplo 1 - Agua de crema de coco
%  min 2x1 + 3x2
%
%  s.a. 0.5x1 + 0.25x2 <= 4
%          x1 +    3x2 >= 20
%          x1 +     x2 == 10
%
%       x1>=0, x2>=0.

c = [ 2  3 ]';
A = [0.5 0.25; 1 3 ; 1 1];
b = [4 20 10]';
inq = [-1 1 0];

p = granM(A,b,c,inq,'min'); p.solve; 

%% Ejemplo 2 - METALCO
%  min 0.4x1 + 0.5x2
%
%  s.a. 0.5x1 + 0.5x2 == 6
%       0.6x1 + 0.4x2 >= 6
%       0.3x1 + 0.1x2 <= 2.7
%
%       x1>=0, x2>=0.

c = [ 0.4  0.5 ]';
A = [0.5 0.5; 0.6 0.4 ; 0.3 0.1];
b = [6 6 2.7]';
inq = [0 1 -1];

p = granM(A,b,c,inq,'min'); p.solve; 

%% Ejemplo 3 - Artesano de Juguetes
%  min 0.4x1 + 0.5x2
%
%  s.a. 0.5x1 + 0.5x2 == 6
%       0.6x1 + 0.4x2 >= 6
%       0.3x1 + 0.1x2 <= 2.7
%
%       x1>=0, x2>=0.

c = [ 1.6  1.4 ]';
A = [10 20; 15 10 ; 18 6];
b = [8000 6000 6300]';
inq = [-1 -1 -1];

p = granM(A,b,c,inq,'max'); p.solve; 

