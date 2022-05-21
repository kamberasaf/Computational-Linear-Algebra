% Task 1

funcs = Functions;
%% section e
num_values = [1,2,5,10,20,50];
h_values = zeros(1,6);
cond_list = zeros(1,6);        % condition number for task a
rel_error_list_b = zeros(1,6); % relative errors for task b
rel_error_list_c = zeros(1,6); % relative errors for task c
rel_error_list_d = zeros(1,6); % relative errors for task d

for k = 1:6
    %% section a
    h=(num_values(k)*pi*Rho)/M;
    h_values(k)=h;
    A = funcs.Matrix_A(M,h);
    V = A*q ;
    [L,U,P] = lu(A);
    K=cond(A);                 % condition number K(A)
    cond_list(k)=cond(A);
    q_norm = norm(q);          % 2-norm of q (p=2).
    v_norm = norm(V);          % 2-norm of V.
    A_fb_norm = norm(A,"fro"); % Frobenius norm of matrix A.
    %% section b
    % we have Aq=V -> LUq=V -> (2) Uq=y -> (1) Ly=V
    y_b = Ly_b(L,P*V);   % calculate y in Ly=V (1)
    q_b = Ux_y(U,y_b);   % calculate q in Uq=y (2)
    rel_q_b_error = funcs.rel_error(q,q_b',2); % calculate relative error.
    rel_error_list_b(k) = rel_q_b_error; % save for graph.

    %% section c
    delta_v=1e-3.*v_norm;
    new_v = V + delta_v;

    y_c = Ly_b(L,P*new_v);   % calculate y in Ly=V (1)
    q_c = Ux_y(U,y_c);       % calculate q in Uq=y (2)
    rel_q_c_error=funcs.rel_error(q,q_c',2); % calculate relative error.
    rel_error_list_c(k) = rel_q_c_error;

    %% section d
    delta_A=A_fb_norm.*1e-3;
    A_d = delta_A + A;
    [L_d,U_d,P_d] = lu(A_d);

    y_d = Ly_b(L_d,P*V);
    q_d = Ux_y(U_d,y_d);
    rel_q_d_error=funcs.rel_error(q,q_d',2);
    rel_error_list_d(k) = rel_q_d_error;
end

% graphs for section e
figure(1);
subplot(2,1,1);
loglog(h_values,rel_error_list_d,LineWidth=1.5); % A
hold on;
loglog(h_values,rel_error_list_c,LineWidth=1.5); % V
loglog(h_values,rel_error_list_b,LineWidth=1.5); % LU
legend('Relative Error _[_A_]','Relative Error _[_V_]','Relative Error _[_L_U_]',Location='northwest')
hold off
xlabel('h_i')
ylabel('Relative Error')
grid on;
title('Relative Error');

subplot(2,1,2);
loglog(h_values,cond_list,LineWidth=1.5,Color='#7E2F8E');
xlabel('h_i')
ylabel('K(A)')
legend('K(A)',Location='northwest')
grid on;
title('Condition Number K(A)');

movegui(1,"southwest");
sgtitle('Question 1: Gaussian Elimination & LU Decomposition')

%% Functions
function y = Ly_b(L,b)  % Solving Lower Triangular System Ly=b.
M = length(b);
y = zeros(1,M);     % result vector
y(1) = b(1)/L(1,1);
sum = 0;            % temp for sigma
for k = 2:M
    sum = b(k);
    for i = 1:k-1
        sum = sum - L(k,i)*y(i);
    end
    y(k)=sum/L(k,k);
    sum=0;
end
end

function x = Ux_y(U,y) % Solving Upper Triangular System Ux=y.
M = length(U);
x = zeros(1,M);    % solution vector
x(M) = y(M)/U(M,M);
sum = 0;           % temp for sigma
for r = M-1:-1:1
    sum = y(r);
    for i = M:-1:r
        sum = sum-U(r,i)*x(i);
    end
    x(r) = sum/U(r,r);
    sum = 0;
end
end
