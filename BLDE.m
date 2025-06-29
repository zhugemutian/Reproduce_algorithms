clc;
clear;
%% 论文相关信息
% Yu Chen, Weicheng Xie, Xiufen Zou. A binary differential evolution algorithm ...
%learning from explored solutions [J]. Neurocomputing, 2015, 149: 1038-1047.

Ztt = 0;
ZZ = [];
ZZ_con_1 = [];
tic
for kkk = 1:30
%% 背包问题参数
load('kp_sc_1000');
% the Instances should include
% profits of all items $P = [p_1,p_2,\dots,p_n]$
% weights of all items $W = [w_1,w_2,\dots,w_n]$
% knapsack capacity $C$
% number of items $n$

E1 = [P./W;1:n];
E2 = sortrows(E1',1,'descend')';
E = E2(2,:);
P = P(E);
W = W(E);
%% BLDE
%算法参数
%最大迭代次数
MaxIter = 200;
%种群个体数量
pop = 30;

pro = 0.5*C/sum(W);
%(Alg. 1, line 1)
X = (rand(pop,n)<pro).*1;

A = (rand(pop,n)<pro).*1;

%变异概率
pm = max(0.05,min(0.15,10/n));

%去除不可行解
for i = 1:pop
    if X(i,:)*W'>C
        X(i,:) = zeros(1,n);
    end
    if A(i,:)*W'>C
        A(i,:) = zeros(1,n);
    end
end

Record = [];

t = 1;
while t <= MaxIter
    %计算X_gb (Alg. 1, line 3)
    T11 = X*P';
    [a,b] = max(T11);
    X_gb = X(b,:);
    %记录
    Record = [Record,a];
    for i = 1:pop
        %选择x,y from X, z from A (Alg. 1, lines 5-6)
        T21 = randperm(pop,2); % x,y
        Rand_x = X(T21(1),:);
        Rand_y = X(T21(2),:);
        T22 = randi(pop); % z
        Rand_z = A(T22(1),:);
        if Rand_y*P'>= Rand_z*P'
            tx = Rand_y;
        else
            tx = Rand_z;
        end
        for j = 1:n
            if Rand_y(j) == Rand_z(j)
                if X_gb(j) ~= Rand_x(j)
                    tx(j) = X_gb(j);
                elseif rand<= pm
                    if rand<= 0.5
                        tx(j) = 1-tx(j); %这里和论文略有区别(Alg. 1, line 13)
                    end
                end
            end
        end
        %判断解tx是否为可行解
        if tx*W' > C %(Alg. 1, lines 18-20)
            k = n;
            while tx*W' > C
                tx(k) = 0;
                k = k-1;
            end
        end
        if tx*P'> X(i,:)*P'
            X(i,:) = tx;
        end
    end
    t = t+1;
    A = X;
end

ZZ = [ZZ,max(Record)];
ZZ_con_1 = [ZZ_con_1;Record];
end
Ztt = toc; %计算时长
ZZ_con = sum(ZZ_con_1/kkk); %收敛曲线
max(ZZ) %多次独立计算的best solution