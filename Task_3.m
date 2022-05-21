%% Task 3

funcs = Functions;

%% task a+b
num_values = [10,5,2,1/2,1/5];
det_A_plot=zeros(1,5);
rel_err=zeros(1,5);
h_values = zeros(1,5);
for i = 1:5
    h=(num_values(i)*pi*Rho)/M;
    h_values(i)=h;
    A = funcs.Matrix_A(M,h);
    det_A_plot(i) = abs(det(A));
    v=A*q;
    A_T = transpose(A);
    pseudo_inv=inv(A_T*A)*A_T;
    rel_q = pseudo_inv * v;
    rel_err(i)=funcs.rel_error(q,rel_q,2);
end

figure(3);
subplot(2,1,1);
loglog(h_values,det_A_plot,'*-',LineWidth=1.5,Color='#77AC30');
title('Q3 Least Squares Solution: Abs(Det(A))');
ylabel('Abs(Det(A))')
xlabel('h')
legend('Relative Error q_a - q',Location='south')
grid on;

subplot(2,1,2);
loglog(h_values,rel_err,'*-',LineWidth=1.5,Color='#D95319');
title('Q3 Least Squares Solution: Relative Error');
xlabel('h')
ylabel('Relative Error')
legend('Relative Error q_a - q',Location='south')
grid on;

movegui(3,"southeast");
sgtitle('Question 3: Least Squares Solutions')