clc;
clear;
%% ���������Ϣ
% Yu Chen, Weicheng Xie, Xiufen Zou. A binary differential evolution algorithm ...
%learning from explored solutions [J]. Neurocomputing, 2015, 149: 1038-1047.

Ztt = 0;
ZZ = [];
ZZ_con_1 = [];
tic
for kkk = 1:30
%% �����������
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
%�㷨����
%����������
MaxIter = 200;
%��Ⱥ��������
pop = 30;

pro = 0.5*C/sum(W);
%(Alg. 1, line 1)
X = (rand(pop,n)<pro).*1;

A = (rand(pop,n)<pro).*1;

%�������
pm = max(0.05,min(0.15,10/n));

%ȥ�������н�
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
    %����X_gb (Alg. 1, line 3)
    T11 = X*P';
    [a,b] = max(T11);
    X_gb = X(b,:);
    %��¼
    Record = [Record,a];
    for i = 1:pop
        %ѡ��x,y from X, z from A (Alg. 1, lines 5-6)
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
                        tx(j) = 1-tx(j); %�����������������(Alg. 1, line 13)
                    end
                end
            end
        end
        %�жϽ�tx�Ƿ�Ϊ���н�
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
Ztt = toc; %����ʱ��
ZZ_con = sum(ZZ_con_1/kkk); %��������
max(ZZ) %��ζ��������best solution