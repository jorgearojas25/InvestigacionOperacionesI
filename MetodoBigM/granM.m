classdef granM < handle
    % granM This class implements big M Simplex Method to solve a
    % linear programming problem in the following format
    % min/max c'x
    % s.t.   Ax {>=, =, <=} b
    % x >= 0
    %
    % Este codigo es una verifiacion de ejercicios hechos a mano entrgados
    % el dia 17 de septiembre del 2018 para la clase de IO, no esta
    % diseñado para resolver probelmas grandes, no es lo suficiente para
    % resolver problemas de informatica de alto rendimiento.
    %
    % Ejcutar example.m
    %
    % 17 Sep 2018
    
    % Jorge Rojas
    % Universidad Distrital Frnacisco Jose de Caldas
    % Basado en el repositorio de:
    % Yiming YAN 
    % University of Edinburgh
    
    properties(SetAccess = private)
        A;                      % Matriz en Ax {<=, =, >=} b
        b;                      % Valores de las restricciones
        c;                      % Vector de coeficientes de la funcion objetivo 
        inq;                    % vector de asignacion, 
                                % 1 para >=; 0 para ==; -1 para <=.
        m;                      % Numero de restricciones
        n;                      % Numero de variables originales
        x_opt;                  % Solucion optima
        fval;                   % Valor optimo de la f objetivo
        type;                   % Minimizar o maximizar 
    end
    
    properties( Access = private )
        M;                     
        
        n_hold;                 % Numero de variables de Holgura añadidas
        n_artfVar;              % Numero de variables artificiales añadidas
        artfcl_var_set          % Indicador de las variables artificiales
        artfcl_var_removed = 1; % 1, if removida
        
        counter = 1;            % Contador de Iteraciones
        
        basis;                  % Vector de indices de variables basicas
        reducedCost;             
        whichRow;                
        
        tab;                    % Tabla
        status;                 % MáxIter, ilimitado, infactible, óptimo
        terminate;
        vars_list;              % Lista de nombres de las variables
        nolimitado_vars;         % Variables ilimitadas
    end
    
    properties( Constant )
        maxIter = 20;
    end
    
    methods
        
        function p = granM(A,b,c,inq, type)
            % CONSTRUCTOR
            
            % Validar entradas
            if nargin < 4
                error('granM: No hay suficientes entradas');
            end
            
            %% Leer la entrada
            % Ver si es max o min
            if strcmpi(type,'max')
                c = -c;
            end
            
            % Validar que b no sea negativa
            b_Neg = b < 0;
            inq(b_Neg) = -inq(b_Neg);
            b(b_Neg) = -b(b_Neg);
            A(b_Neg,:) = -A(b_Neg,:);
            
            % Asignacion de propiedades
            p.A = A; p.b = b; p.c = c;
            p.inq = inq;
            p.terminate = 0;
            
            % Asignar el valor de M
            p.M = ceil(max(norm(A),norm(b))/100)*100;
            
            [p.m, p.n] = size(A);
            
            % Crear la lista de variables
            p.vars_list = cell(1,p.n);
            for i=1:p.n
                p.vars_list{i} = ['x' num2str(i)];
            end
            
            p.type = type;
        end
        
        %% Funcion principal
        function solve(p)
            % Esta es la funcion que resuelve el problema.
            
            p.transform_to_standardForm;
            p.construct_initial_tableau;
            p.printTab;
            while ~p.terminate
                p.do_simplex_iterate;
                p.increment_counter;
                p.printTab;
                p.remove_artfcl_vars;
                p.check_termination;
            end
            p.output_summary;
        end
        
        function transform_to_standardForm(p)
            
            % Añade las variables de holguta
            p.n_hold = sum(p.inq < 0)+sum(p.inq > 0);
            p.A = [p.A zeros(p.m, p.n_hold)];
            p.c = [p.c; zeros(p.n_hold,1)];
            
            idx_slack = p.n+1;
            for i = 1:p.m
                switch p.inq(i)
                    case 1  % >=
                        p.A(i, idx_slack) = -1;
                        idx_slack = idx_slack+1;
                         p.vars_list{end+1} = ['s' num2str(i)];
                    case -1 % <=
                        p.A(i, idx_slack) = 1;
                        idx_slack = idx_slack+1;
                        p.vars_list{end+1} = ['s' num2str(i)];
                end
            end
        end
        
        function construct_initial_tableau(p)

            % Encuentra la base potencial
            for i = 1:p.n+p.n_hold
                
                if nnz(p.A(:,i)) == 1
                    row_number = find(p.A(:,i) == 1);
                    if ~isempty(row_number)                        
                        % Revisa que la fila actual sea evaluada
                        if ~ismember(row_number, p.whichRow)
                            p.whichRow(end+1) = row_number;
                            p.basis(row_number) = i;
                        end
                    end
                end
            end
            
            % Añade las M's necesarias            
            p.n_artfVar = p.m - sum(p.basis > 0);
            if  p.n_artfVar > 0
                p.artfcl_var_removed = 0;
                
                % Registra el indice de variables artificiales
                p.artfcl_var_set = (p.n + p.n_hold + 1) : (p.n + p.n_hold + p.n_artfVar);
                
                p.A = [p.A zeros(p.m,p.n_artfVar)];
                
                add_to_rows = setdiff(1:p.m, p.whichRow);
                
                % Formula la matriz A con las variables artificiales
                % agregadas
                for i = 1:length(add_to_rows)
                   p.A(add_to_rows(i), p.n+p.n_hold+i) = 1;
                   p.vars_list{end+1} = ['a' num2str(i)];
                   p.basis(add_to_rows(i)) = p.n+p.n_hold+i;
                end
                
                p.c = [p.c; p.M*ones(p.n_artfVar,1)];
                
            end
                        
            % Recibe las bases
            p.reducedCost = zeros(1,p.n + p.n_hold + p.n_artfVar);
            y = p.c(p.basis);
            
            p.reducedCost = (p.c - p.A'*y)';
            
            % Construye la tabla inicial
            p.tab = [p.A  p.b];
            p.tab = [p.tab; [ p.reducedCost -p.b'*y ]];
        end
        
        function do_simplex_iterate(p)
            % DO_SIMPLEX_ITERATE Perform one iterate of simplex method.
            
            % Encuentra el valor mas pequeño
            [~, workCol] = min( p.tab(end,1:end-1) );
            
            % elige el residuo mas pequeño
            work_col_positive_elems = find( p.tab(1:end-1,workCol) > 0 );
            
            [~, chooseRow] = min( p.tab(work_col_positive_elems, end) ./...
                p.tab(work_col_positive_elems, workCol) );
            chooseRow = work_col_positive_elems(chooseRow);
            
            % Entrada y salida de las bases
            p.basis(chooseRow) = workCol;
            
            % Recibe el pivot
            pivot = p.tab(chooseRow, workCol);
            
            % Elimina elementos trabajados
            % y añade el elemento actual a ala tabla.
            p.tab(chooseRow, :) = p.tab(chooseRow,:) / pivot;
            
            tmp_row_indx = setdiff(1:p.m+1,chooseRow);
            p.tab(tmp_row_indx, :) = p.tab(tmp_row_indx, :) -...
                ones(p.m,1) * p.tab(chooseRow,:) .*...
                ( p.tab(tmp_row_indx, workCol) *...
                ones(1,size(p.tab,2)) );
        end
        
        function remove_artfcl_vars(p)
            % Mira que las variables artificiales sean removidas
            if ~p.artfcl_var_removed
                if all( p.tab(end, p.artfcl_var_set) > 0 )
                    p.tab(:,p.artfcl_var_set) = [];
                    p.vars_list( p.artfcl_var_set ) = [];
                    p.artfcl_var_removed = 1;
                    p.printTab;
                end
            end
        end
        
        function check_termination(p)
            p.terminate = 0;
            
            % Evita que supere el maximo de iteraciones
            if p.counter >= p.maxIter
                p.terminate = 1;
                p.status = 'maxIter';
            end
            
            % Revisa los productos de M que no sean negativos
            if all(p.tab(end,1:end-1) > -1e-10)
                p.terminate = 1;
                
                % No factible
                if ~p.artfcl_var_removed
                    p.status = 'infactible';
                else
                    % Optimo
                    p.status = 'optimo';
                    p.x_opt = zeros(p.n + p.n_hold + p.n_artfVar, 1);
                    p.x_opt(p.basis) = p.tab(1:end-1,end);
                    if strcmpi(p.type,'max')
                        p.fval = p.tab(end,end);
                    elseif strcmpi(p.type,'min')
                        p.fval = -p.tab(end,end);
                    end
                end
            else
                % Ilimitado
                negative_rc_cols = find(p.tab (end,1:end-1) < 0);
                nolimitado_idx = sum(p.tab(1:end-1, negative_rc_cols) < 0);
                nolimitado_idx = nolimitado_idx == p.m;
                if any(nolimitado_idx)
                    p.status = 'nolimitado';
                    p.terminate = 1;
                    p.nolimitado_vars = negative_rc_cols(nolimitado_idx > 0);
                end                
            end
        end
        
        
        function increment_counter(p)
            p.counter = p.counter+1;
        end
        
        function printTab(p)
            fprintf('\n\n');
            fprintf('\t\tTabla %2s:\n\n', num2str(p.counter));
            % Print header
            for k=1:length(p.vars_list)+1
                if k==1
                    fprintf('%2s | ', '');
                else
                    fprintf('%8s %4s', p.vars_list{k-1}, '');
                    if k == length(p.vars_list)+1
                        fprintf(' | %8s', 'RHS');
                    end
                end
            end
            fprintf('\n');
            fprintf('------------------------------------------------------------------\n');
            
            for i=1:p.m+1
                if i <= p.m
                    fprintf('%2s | ', p.vars_list{p.basis(i)});
                else
                    fprintf('------------------------------------------------------------------\n');
                    fprintf('%2s | ', 'RC');
                end
                for j = 1:size(p.tab,2)
                    fprintf('%8.2f %4s', p.tab(i,j), '');
                    if j == size(p.tab,2)-1
                       fprintf(' | ') 
                    end
                end
                fprintf('\n');
            end
            
        end
        
        function output_summary(p)
            if p.terminate
                fprintf('\n========== ========== ==========\n');
                fprintf('Terminado. \n');
                switch lower( p.status )
                    case 'optimo'
                        fprintf('Soluciones optimas: \n');
                        for i=1:p.n
                            fprintf('x[%2d ] = %5.2f\n', i, p.x_opt(i));
                        end
                        fprintf('Valor Optimo de la funcion objetivo: %5.2f\n', p.fval);
                        
                    case 'nolimitado'
                        fprintf('El problema no esta limitado. \n');
                        fprintf('Variables no limitadas: ');
                        
                        for i = 1:length(p.nolimitado_vars)
                            fprintf('%2s ',p.vars_list{p.nolimitado_vars(i)});
                        end
                        fprintf('\n');
                        
                    case 'infactible'
                        fprintf('El problema es infactible.\n');
                    case 'maxiter'
                        fprintf('Alcanzo el maximo numero de iteraciones posibles.\n');
                end
                
            end
        end
    end
end