classdef Functions
    % This file is for functions used in all 3 Tasks.
    methods
        function res = Matrix_A(~,M,h)
            % create matrix A with given M (size) and h (changes values).
            res = zeros(M);
            Rho = 1;
            for m = 1:M
                for n = 1:M
                    rad = sqrt((h+Rho*sin((m*pi)/M)-Rho*sin((n*pi)/M)).^2+ ...
                        (Rho*cos((m*pi)/M)-Rho*cos((n*pi)/M)).^2);
                    res(m,n) = 1 ./ (4*pi*rad);
                end
            end
        end

        function res = Matrix_A_task2(~,M,h) % create matrix A with given M (size) and h (changes values).
            res = zeros(M);
            Rho = 1;
            for m = 1:M
                for n=1:M
                    rad = (h+Rho*sin((m*pi)/M)-Rho*sin((n*pi)/M)).^2+(Rho*cos((m*pi)/M)-Rho*cos((n*pi)/M)).^2;
                    res(m,n)=1 ./ (4*pi*rad);
                end
            end
        end

        function res = abs_error(~,e,e_t,p) % get absolute error.
            res = norm(e-e_t,p);
        end

        function [res,per] = rel_error(~,e,e_t,p) % get relative error.
            res = norm(e-e_t,p)/norm(e,p); % relative error.
            per = 100 * res; % percent error.
        end

        function [iter_track,dist_error,real_error] = Gauss_Seidel(~,A,x,b)
            tol = 1e-3;
            L = tril(A,-1);    % Lower triangular matrix
            U = triu(A,1);     % Upper triangular matrix
            D = diag(diag(A)); % Diagonal matrix
            Q = L+D;
            G = -inv(Q)*U;
            c = inv(Q)*b;
            itr = 1; max_iter = 600;
            iter_track = ones;
            q_k = c;           % q^(0) = 0 , q^(1) = c
            dist_error = ones; % ||q^(1)-q^(0)|| / ||q^(0)|| (can't divide by 0)
            real_error = ones; % ||q^(1)-q|| / ||q||
            real_error(1) = norm(c-x,2)/norm(x);
            error = norm(x-q_k(:, itr),"inf")/norm(x,"inf");
            while abs(error)>tol && itr<=max_iter        % make sure we stop.
                q_k(:,itr+1) = G*q_k(:,itr) + c;                 % Gauss-Seidel Algorithm.
                error = norm(x-q_k(:, itr),"inf")/norm(x,"inf"); % finding error ||q^k-q||
                dist_error(itr+1) = norm(q_k(:,itr+1) - q_k(:, itr), ...
                    "inf")/norm(q_k(:, itr),"inf");
                real_error(itr+1) = norm(q_k(:,itr+1) - x, ...
                    "inf")/norm(x,"inf");
                iter_track(itr+1)=itr+1;
                itr = itr + 1;
            end
        end


        function [iter_track,dist_error,real_error] = Jacobi(~,A,x,b) %Aq=v,Ax=b
            tol = 1e-3;
            D = diag(diag(A));        % Diagonal matrix
            Q = D;
            I = eye(length(x));
            G = I-inv(Q)*A;
            norm_G = norm(G,"inf");   % Check for || G || < 1,
            c = inv(Q)*b;
            itr = 1; max_iter = 600;
            iter_track = ones;
            q_k = c;           % q^(0)=0 , q^(1)=c
            dist_error = ones; real_error = ones;  % same as G-S
            real_error(1) = norm(c-x,2)/norm(x);   % ".
            error = norm(x-q_k(:, itr),"inf")/norm(x,"inf"); % just a start value.
            while abs(error) > tol && itr<=max_iter % make sure we stop.
                q_k(:,itr+1) = G*q_k(:,itr) + c;         % Jacobi Algorit.
                error = norm(x-q_k(:, itr),"inf")/norm(x,"inf");   % finding error
                dist_error(itr+1) = norm(q_k(:,itr+1) - q_k(:, itr), ...
                    "inf")/norm(q_k(:, itr),"inf");
                real_error(itr+1) = norm(q_k(:,itr) - x,"inf")/norm(x,"inf");
                iter_track(itr+1)=itr+1;
                itr = itr + 1;
            end
        end

    end
end