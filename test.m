%% ABABIO GODFRED OPOKU AND VARUN 
[elk] = yourinitials_estiff (10, 100, 20, 1, 8, 2, 29000, 0.3, 120*sqrt(2));

[gamma] = yourinitials_etran([0,0,0], [120,120,0], [-0.7071,0.7071,0]);

K = gamma'*elk*gamma;