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

%% Example 3 - check adding big M

c = [1 2 3]';
A = [0 1 2; 3 2 1];
b = [3 6]'; 
inq = [0 0 0];
p = granM(A,b,c,inq,'min'); p.solve; 

%% Example 4 - check infeasible
%  max   x1 + 2x2
% s.t.   x1 +  x2 >= 4
%       3x1 + 2x2 <= 6
%        x1, x2 >=0
c = [1 2]';
A = [1 1 ; 3 2];
b = [4 6]';
inq = [1 -1];

p = granM(A,b,c,inq,'max'); p.solve; 

%% Example 5 - check unbounded
%  min -x1-2x2
% s.t.    x1 + x2 >= 3
%         x1 - x2 >= 1
%       x1, x2 >= 0

c = [-1 -2]';
A = [1 1; 1 -1];
b = [3 ; 1];
inq = [1 1];

p = granM(A,b,c,inq,'min'); p.solve; 

%% Example 6 - check infeasible with granM

c = [ 1 1 1 1 ]';
A = [ 1 2 0 1; 0 0 1 3; 1 0 1 0 ];
b = [2 3 6 ]';
inq = [ -1 -1 1 ];
p = granM(A,b,c,inq,'max'); p.solve; 

%% Example 7 - check unbounded with granM
c = [ 1 2 3 4 ]';
A = [ 1 2 -1 2; 2 0 3 -1 ];
b = [ 4 6 ]';
inq = [ 1 0 ];
p = granM(A,b,c,inq,'max'); p.solve; 




