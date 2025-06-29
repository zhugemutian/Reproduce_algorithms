clc;
clear;
%% 论文相关信息
% Hu Peng Zhijian Wu, Peng Shao, Changshou Deng. Dichotomous Binary Differential...
%Evolution for Knapsack Problems [J]. Mathematical Problems in Engineering, 2016, 5732489.

%1807 3403 5444 9495 18844
% 659 1332 1963 3250 6482
% 813 1631 2433 4078 8228
Ztt = 0;
ZZ = [];
ZZ_con_1 = [];
tic
for kkk = 1:30
%% 背包问题参数
load('kp_sc_1000');
%the Instances should include
% profits of all items $P = [p_1,p_2,\dots,p_n]$
% weights of all items $W = [w_1,w_2,\dots,w_n]$
% knapsack capacity $C$
% number of items $n$
E1 = [P./W;1:n];
E2 = sortrows(E1',1,'descend')';
E = E2(2,:);
P = P(E);
W = W(E);

%% DBDE
%算法参数
%最大迭代次数
MaxIter = 200;
%种群个体数量
pop = 30;

%CR
cr1 = 0.2;
cr2 = 0.5;

%% 算法开始
%X = randi([0,1],pop,n); %(lines 1-7)
X = (rand(pop,n)>0.2)*1;
for i = 1:pop
    j = n;
    while X(i,:)*W'>C
        X(i,j) = 0;
        j = j-1;
    end
end

Record = [];

t = 1;
while t <= MaxIter
    for i = 1:pop
        %Dichotomous Mutation
        T1 = randperm(pop,2);
        t11 = T1(1);
        t12 = T1(2);
        for j = 1:n
            if X(t11,j) == X(t12,j)
                Vig(j) = X(t11,j);
            else
                if rand<0.5
                    Vig(j) = 1;
                else
                    Vig(j) = 0;
                end
            end
            %Dichotomous Crossover
            %决定CR数值
            if X(t11,j) == X(t12,j)
                cr = cr1;
            else
                cr = cr2;
            end
            if rand <= cr
                Uig(j) = Vig(j);
            else
                Uig(j) = X(i,j);
            end
        end
        for i = 1:pop
            j = n;
            while Uig*W'>C
                Uig(j) = 0;
                j = j-1;
            end
        end
        if Uig*P'>X(i,:)*P'
            X(i,:) = Uig;
        end
    end
    [a,b] = max(X*P');
    Record = [Record, a];
    t = t+1;
end
%plot(1:MaxIter,Record)

ZZ = [ZZ,max(Record)];
ZZ_con_1 = [ZZ_con_1;Record];
end
Ztt = toc;
ZZ_con = sum(ZZ_con_1/kkk);
max(ZZ)