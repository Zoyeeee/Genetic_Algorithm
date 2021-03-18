clear all;clc;
NP=100;%种群数量
NG=500;%代数
Pc=0.8;%杂交概率
Pm=0.15;%变异概率
NZ=30;%种群内一个样本，样本中含个体的数量
q=0.3;%最好样本的选择概率
%%
%初始化第一代
X=zeros(NZ,NP);
Y=zeros(1,NP);
plot_Y=zeros(1,NG);
plot_X=1:NG;
for j=1:NP
    for i=1:NZ
        X(i,j)=rand()*10-5;
    end
end

for g=1:NG
    for j=1:NP
        Y(j)=fitness(X(:,j),NZ);
    end
    %%
    %选择一个好样本
    [Y,id_aftersort]=sort(Y,'descend');
    X=X(:,id_aftersort);
    plot_Y(g)=Y(1);%大循环为g=1:NG
    Y(1)=q;

    for i=2:NP
        Y(i)=(q*((1-q)^(i-1)))/(1-((1-q)^NP));%各样本选择的概率
    end
    for i=2:NP
        Y(i)=Y(i-1)+Y(i);
    end

    F_X=X;
    for k=1:NP
        sita=rand();
        for i=1:NP
            if sita<=Y(i)
                Select_sample_father=i;
                break
            end
        end
        Select_sample_mather=floor(rand()*(NP-1))+1;

        for i=1:NZ
            Putcross=rand();
            if Putcross<=Pc
                alpha=rand()*1.5-0.25;
                X(i,k)=F_X(i,Select_sample_father)+alpha*(F_X(i,Select_sample_mather)-F_X(i,Select_sample_father));
            end
            Putvariation=rand();
            if Putvariation<=Pm
                %高斯变异
                mean_i=mean(F_X(i,:));
                cov_i=cov(F_X(i,:));
                X(i,k)=mvnrnd(mean_i,cov_i,1);
            end
        end
    end
end

%%
%画图
plot(plot_X,plot_Y,'g','LineWidth',1);
hold on;
plot(plot_X,plot_Y,':r','LineWidth',3);
title('每一代最大值曲线');
xlabel('代数');
ylabel('最大值')
hold off;
%%
%适应度函数基本定义
%基本数据定义
function y=fitness(X,NZ)
pi=3.14;
C1=20;
C2=0.2;
C3=2*pi;
e=2.71282; %如果和精度相关得注意，有可能有效位数超过了
cusum1=X'*X;
cusum2=0;
for i=1:NZ
    cusum2=cos(C3*X(i))+cusum2;
end
y=C1*exp(-C2*sqrt(0.5*cusum1))+exp(0.5*cusum2)-C1-e;
end
%%
