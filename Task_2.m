%% Task 2

funcs = Functions;
Rho=1;

iter_list_5M = zeros;
dist_error_5M = zeros;
rel_error_5M = zeros;
iter_list_2M = zeros;
dist_error_2M = zeros;
rel_error_2M = zeros;
iter_list_1M = zeros;
dist_error_1M = zeros;
rel_error_1M = zeros;

%% task a + b
num_values_2 = [5,2,1];
for i = 1:3
    h = (pi*Rho)/(num_values_2(i)*M);
    A = funcs.Matrix_A(M,h);
    v = A*q;
    [iter_list,dist_error,rel_error] = funcs.Gauss_Seidel(A,q,v);

    if i==1
        iter_list_5M = iter_list;
        dist_error_5M = dist_error;
        rel_error_5M = rel_error;
    end

    if i==2
        iter_list_2M = iter_list;
        dist_error_2M = dist_error;
        rel_error_2M = rel_error;
    end
    if i ==3
        iter_list_1M = iter_list;
        dist_error_1M = dist_error;
        rel_error_1M = rel_error;
    end

end

%% task c
h=(pi*Rho)/(5*M);
A=funcs.Matrix_A(M,h);
v=A*q;
[iter_list_c,dist_error_c,rel_error_c] = funcs.Jacobi(A,q,v);



%% task d
A_d = funcs.Matrix_A_task2(M,h);
v=A_d*q;
[iter_list_d,dist_error_d,rel_error_d] = funcs.Jacobi(A_d,q,v);

Q2 = figure('Renderer', 'painters', 'Position', [13 11 700 550]);
hold on
subplot(13,11,[1,38]) % Gauss-Seidel 5M
semilogy(iter_list_5M,dist_error_5M,iter_list_5M,rel_error_5M);
legend('q^(^k^)-q^(^k^-^1^)','q^(^k^)-q',Location='southwest');
title('G-S Relative Errors (5M)');
ylabel('Relative Errors');
xlabel('Iterations');
xlim([1,iter_list_5M(end)])

subplot(13,11,[7,44]) % Gauss-Seidel 2M
semilogy(iter_list_2M,dist_error_2M,iter_list_2M,rel_error_2M);
legend('q^(^k^)-q^(^k^-^1^)','q^(^k^)-q',Location='southwest');
ylabel('Relative Errors');
xlabel('Iterations');
title('G-S Relative Errors (2M)');
xlim([1,iter_list_2M(end)])

subplot(13,11,[59,85]) % Gauss-Seidel 1M
semilogy(iter_list_1M,rel_error_1M,iter_list_1M,dist_error_1M);
ylabel('Relative Errors');
xlabel('Iterations');
legend('q^(^k^)-q','q^(^k^)-q^(^k^-^1^)',Location='east');
title('G-S Relative Errors (1M)')


subplot(13,11,[100,137]) % Jacobi (c)
semilogy(iter_list_c,rel_error_c,iter_list_c,dist_error_c);
ylabel('Relative Errors');
xlabel('Iterations');
legend('q^(^k^)-q','q^(^k^)-q^(^k^-^1^)',Location='northwest');
title('Jacobi Relative Errors (5M)_C')


subplot(13,11,[106,143]) % Jacobi (d)
semilogy(iter_list_d,dist_error_d,iter_list_d,rel_error_d);
legend('q^(^k^)-q^(^k^-^1^)','q^(^k^)-q',Location='northeast');
title('Jacobi Relative Errors (5M)_D');
ylabel('Relative Errors');
xlabel('Iterations');

sgtitle('Question 2: Iteration Solutions')
movegui(Q2,"north");