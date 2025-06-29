clc;
clear;
% Kennedy, J.; Eberhart, R.C. Discrete binary version of the particle swarm...
% algorithm. In Proceedings of the IEEE International Conference on Systems, ...
% Man and Cybernetics, 1997, 5: 4104–4108.

% 论文部分地方根据下面的文献进行了适当的修改
% Hu Peng Zhijian Wu, Peng Shao, Changshou Deng. Dichotomous Binary Differential...
%Evolution for Knapsack Problems [J]. Mathematical Problems in Engineering, 2016, 5732489.

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
%% BPSO 
%算法参数
%最大迭代次数
MaxIter = 200;
%学习因子
c1 = 2;
c2 = 2;
%种群个体数量
pop = 30;
%惯性系数
%omega = 0.5;%也可以随着迭代次数逐渐减小
%算法初始化
%速度上界
Vmax = 6; %避免速度过大或过小导致陷入局部最优解无法跳出
%速度初始化
%速度分布为[-Vmax, Vmax]
V = rand(pop,n)*Vmax-1.6*Vmax; %多减一点, 减少被选择的物品, 避免出现大多数解为不可行解的情况...
%1.4数值的设定与C和sum(W)的比值相关
%选择概率
SV = 1./(1+exp(-V));
%速度所对应的位置
X = [rand(pop,n)<SV]*1;

for i1 = 1:pop
    if X(i1,:)*W' <= C %剔除不可行解
        pBest_value(i1) = X(i1,:)*P';
        pBest(i1,:) = X(i1,:);
    else %对于不可行解用全0解替换
        pBest_value(i1) = 0;
        pBest(i1,:) = -1*ones(1,n);
    end
end
[a,b] = max(pBest_value);
gBest_value = a;
gBest = X(b,:);
%记录算法迭代
Record = gBest_value;
%算法循环
for Iter = 2:MaxIter
    %速度更新
    for i1 = 1:pop %不一次性迭代是为了每个元素的随机数不同
        for i2 = 1:n
            V(i1,i2) = V(i1,i2)+c1*rand*(X(i1,i2)-pBest(i1,i2))+c2*rand*(X(i1,i2)-gBest(i2));
            %对于速度进行边界修复
            if V(i1,i2)>Vmax
                V(i1,i2) = Vmax;
            end
            if V(i1,i2)<-Vmax
                V(i1,i2) = -Vmax;
            end
        end
    end
	%更新选择概率
    SV = 1./(1+exp(-V));
    %更新位置
    X = [rand(pop,n)<SV]*1;
    for i = 1:pop
        j = n;
        while X(i,:)*W'>C
            if X(i,j) == 1;
                X(i,j) = 0;
            end
            j = j-1;
        end
    end
    %更新个体最优和全局最优
    for i1 = 1:pop
        if X(i1,:)*W' <= C %剔除不可行解
            if X(i1,:)*P' > pBest_value(i1)
                pBest_value(i1) = X(i1,:)*P';
                pBest(i1,:) = X(i1,:);
            else %对于不可行解用全0解替换
                pBest_value(i1) = 0;
                pBest(i1,:) = -1*ones(1,n);
            end
        end
    end
    [a,b] = max(pBest_value);
    gBest_value = max(Record(Iter-1),a);
    gBest = X(b,:);
    Record = [Record,gBest_value];
end
%plot(1:MaxIter,Record)

ZZ = [ZZ,max(Record)];
ZZ_con_1 = [ZZ_con_1;Record];
end
Ztt = toc; %计算时长
ZZ_con = sum(ZZ_con_1/kkk); %收敛曲线
max(ZZ) %多次独立计算的best solution