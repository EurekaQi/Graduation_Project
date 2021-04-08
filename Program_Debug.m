%%
clc,clear
cd 'D:\学习用文件\课题组\项目\毕业设计\程序&数据'
%%
%读取数据
[~,~,BasicPara]=xlsread('IEEE30.xlsm','BasicPara');
[~,~,NetBuses]=xlsread('IEEE30.xlsm','NetBuses');
[~,~,Loads]=xlsread('IEEE30.xlsm','Loads');
[~,~,NetBranches]=xlsread('IEEE30.xlsm','NetBranches');
[~,~,LoadCurves]=xlsread('IEEE30.xlsm','LoadCurves');
[~,~,UnitThermalGenerators]=xlsread('IEEE30.xlsm','UnitThermalGenerators');
[~,~,WindFarmGenerators]=xlsread('IEEE30.xlsm','WindFarmGenerators');
[~,~,PhotovoltaicGenerators]=xlsread('IEEE30.xlsm','PhotovoltaicGenerators');
[~,~,Scenarios]=xlsread('IEEE30.xlsm','Scenarios');
[~,~,BidData]=xlsread('IEEE30.xlsm','BidData');
disp('数据导入完毕');
%%
%基本参数处理
T=BasicPara{2,1};                             %设置时段数
T_Gap=BasicPara{2,2};                         %设置每段时间间隔
W=BasicPara{2,17};                            %设置场景数
Seg_offer=BasicPara{2,3};                     %发电机报价段数
Seg_bid=BasicPara{2,4};                       %用户报价段数
Offer_UpperCup=BasicPara{2,7};                %用户报价上限
Offer_LowerCup=BasicPara{2,8};                %用户报价下限
Gama_FT=BasicPara{2,12};                      %未来TGC价格因数
Theta_FT=BasicPara{2,13};                     %TGC充足时价格因数
Num_Approximate=BasicPara{2,14};              %近似基准值数量
RPS_ratio=BasicPara{2,15};                    %RPS需求比率
Penalty_Price=BasicPara{2,16};                %TGC不足的惩罚价格
RT_Wave=[0.992;1.025;1.012;0.994;1.001;1.017];%实时市场中的波动
RE_Gens=[7;8;10;12;14];                       %新能源发电公司所拥有的机组编号
RE_Gens_Other=[9;11;13];                      %其他新能源发电机编号
M=400000;
gap=1e-5;
Q_max=20;  Q_min=0;
Q_ref_gap=(Q_max-Q_min)/(Num_Approximate-1);                
%节点数据整理
NetBuses_Reform=zeros(size(NetBuses,1)-1,1);
for i=1:size(NetBuses,1)-1
    if length(NetBuses{i+1,2})==2
        NetBuses_Reform(i)=str2double(NetBuses{i+1,2}(2));
    else
        NetBuses_Reform(i)=str2double(NetBuses{i+1,2}(2:3));
    end
end

NetBranches_Reform=zeros(size(NetBranches,1)-1,4);
for i=1:size(NetBranches,1)-1
    if length(NetBranches{i+1,2})==5
        NetBranches_Reform(i,1)=str2double(NetBranches{i+1,2}(5));
    else
        NetBranches_Reform(i,1)=str2double(NetBranches{i+1,2}(5:6)); 
    end
    
    if length(NetBranches{i+1,3})==5
        NetBranches_Reform(i,2)=str2double(NetBranches{i+1,3}(5));
    else
        NetBranches_Reform(i,2)=str2double(NetBranches{i+1,3}(5:6)); 
    end
    
    NetBranches_Reform(i,3)=NetBranches{i+1,4};
    NetBranches_Reform(i,4)=NetBranches{i+1,5};
end
%节点数据处理
AllNodes=cell(size(NetBuses_Reform,1),1);
for i=1:size(NetBuses_Reform,1)
    tmp=[];
    for j=1:length(NetBranches_Reform)   
        if NetBranches_Reform(j,1)==NetBuses_Reform(i)
            tmp=[tmp;NetBranches_Reform(j,2)];
        end
        if NetBranches_Reform(j,2)==NetBuses_Reform(i)
            tmp=[tmp;NetBranches_Reform(j,1)]; %#ok<*AGROW>
        end
    end
    AllNodes{i}=struct('connectnodes',tmp);
end
%线路数据处理
B=zeros(size(NetBuses_Reform,1),size(NetBuses_Reform,1));
P_Lmax=zeros(size(NetBuses_Reform,1),size(NetBuses_Reform,1));
for i=1:size(NetBranches_Reform,1)
    B(NetBranches_Reform(i,1),NetBranches_Reform(i,2))=1/NetBranches_Reform(i,3);
    B(NetBranches_Reform(i,2),NetBranches_Reform(i,1))=1/NetBranches_Reform(i,3);
    P_Lmax(NetBranches_Reform(i,1),NetBranches_Reform(i,2))=NetBranches_Reform(i,4);
    P_Lmax(NetBranches_Reform(i,2),NetBranches_Reform(i,1))=NetBranches_Reform(i,4);
end
%负荷数据处理
AllLoads=cell(size(Loads,1)-1,1);
for i=1:size(Loads,1)-1
    powerconsume_DA=zeros(T,1);
    powerconsume_RT=zeros(T,1);
    for j=1:size(LoadCurves,1)-1
        if rem(j,T_Gap)==0
            powerconsume_DA(j/T_Gap)=Loads{i+1,4}*(LoadCurves{j-2,2}+LoadCurves{j-1,2}+LoadCurves{j,2}+LoadCurves{j+1,2})/4;
            powerconsume_RT(j/T_Gap)=powerconsume_DA(j/T_Gap)*RT_Wave(j/T_Gap);
        end
    end
    
    
    if length(Loads{i+1,2})==5
        location=str2double(Loads{i+1,2}(5));
    else
        location=str2double(Loads{i+1,2}(5:6));
    end
    
    AllLoads{i}=struct('location',location,'powerconsume_DA',powerconsume_DA,'powerconsume_RT',powerconsume_RT); 
end
%火电机组数据处理
AllThermalGenerators=cell(size(UnitThermalGenerators,1)-1,1);
for i=1:size(UnitThermalGenerators,1)-1
    if length(UnitThermalGenerators{i+1,3})==5
        location=str2double(UnitThermalGenerators{i+1,3}(5));
    else
        location=str2double(UnitThermalGenerators{i+1,3}(5:6));
    end
    
    capacity=UnitThermalGenerators{i+1,4};
    G_min=UnitThermalGenerators{i+1,5};
    marign_cost=UnitThermalGenerators{i+1,12};
    
    AllThermalGenerators{i}=struct('location',location,'capacity',capacity,'P_min',G_min,'marign_cost',marign_cost);
end
%风力发电机组数据处理
AllWindFarmGenerators=cell(size(WindFarmGenerators,1)-1,1);
for i=1:size(WindFarmGenerators,1)-1
    if length(WindFarmGenerators{i+1,3})==5
        location=str2double(WindFarmGenerators{i+1,3}(5));
    else
        location=str2double(WindFarmGenerators{i+1,3}(5:6));
    end

    capacity=WindFarmGenerators{i+1,4};
    G_min=WindFarmGenerators{i+1,5};
    generationcost=WindFarmGenerators{i+1,6};

    typicaloutput=zeros(T,1);
    for j=1:T
        typicaloutput(j)=WindFarmGenerators{i+1,j+6}*WindFarmGenerators{i+1,4};
    end

    AllWindFarmGenerators{i}=struct('location',location,'capacity',capacity,'P_min',G_min,'generationcost',generationcost,'typicaloutput',typicaloutput);
end
%光伏发电机组数据处理
AllPhotovoltaicGenerators=cell(size(PhotovoltaicGenerators,1)-1,1);
for i=1:size(PhotovoltaicGenerators,1)-1
    if length(PhotovoltaicGenerators{i+1,3})==5
        location=str2double(PhotovoltaicGenerators{i+1,3}(5));
    else
        location=str2double(PhotovoltaicGenerators{i+1,3}(5:6));
    end

    capacity=PhotovoltaicGenerators{i+1,4};
    G_min=PhotovoltaicGenerators{i+1,5};
    generationcost=PhotovoltaicGenerators{i+1,6};

    typicaloutput=zeros(T,1);
    for j=1:T
        typicaloutput(j)=PhotovoltaicGenerators{i+1,j+6}*PhotovoltaicGenerators{i+1,4};
    end

    AllPhotovoltaicGenerators{i}=struct('location',location,'capacity',capacity,'P_min',G_min,'generationcost',generationcost,'typicaloutput',typicaloutput);
end
%将所有机组数据整合至一个元胞之中
AllGens=[AllThermalGenerators;AllWindFarmGenerators;AllPhotovoltaicGenerators];   
%场景数据处理
AllScenarios=cell(size(Scenarios,1)-1,1);
for i=1:size(Scenarios,1)-1
    probability=Scenarios{i+1,3};
    WF_change=Scenarios{i+1,4};
    PV_change=Scenarios{i+1,5};

    AllScenarios{i}=struct('probability',probability,'WF_change',WF_change,'PV_change',PV_change);
end
%发电机组出力边际成本确定
Marign_Cost=zeros(size(AllGens,1),Seg_offer);
for i=1:size(AllThermalGenerators,1)+size(AllWindFarmGenerators,1)+size(AllPhotovoltaicGenerators,1)
    for j=1:Seg_offer
        if i<=size(AllThermalGenerators,1)
            Marign_Cost(i,j)=AllThermalGenerators{i}.marign_cost;
        elseif i>size(AllThermalGenerators,1)&&i<=size(AllThermalGenerators,1)+size(AllWindFarmGenerators,1)
            Marign_Cost(i,j)=AllWindFarmGenerators{i-size(AllThermalGenerators,1)}.generationcost;
        else
            Marign_Cost(i,j)=AllPhotovoltaicGenerators{i-size(AllThermalGenerators,1)-size(AllWindFarmGenerators,1)}.generationcost;
        end
    end
end
%日前市场发电机每报价段出力上下限确定
P_DA_Gmin=zeros(size(AllGens,1),Seg_offer,T);
P_DA_Gmax=zeros(size(AllGens,1),Seg_offer,T);
for i=1:size(AllThermalGenerators,1)+size(AllWindFarmGenerators,1)+size(AllPhotovoltaicGenerators,1)
    for j=1:Seg_offer
        for k=1:T
            if i<=size(AllThermalGenerators,1)
                P_DA_Gmax(i,j,k)=AllThermalGenerators{i}.capacity;
                P_DA_Gmin(i,j,k)=0;
            elseif i>size(AllThermalGenerators,1)&&i<=size(AllThermalGenerators,1)+size(AllWindFarmGenerators,1)
                P_DA_Gmax(i,j,k)=AllWindFarmGenerators{i-size(AllThermalGenerators,1)}.typicaloutput(k);
                P_DA_Gmin(i,j,k)=0;
            else
                P_DA_Gmax(i,j,k)=AllPhotovoltaicGenerators{i-size(AllThermalGenerators,1)-size(AllWindFarmGenerators,1)}.typicaloutput(k);
                P_DA_Gmin(i,j,k)=0;
            end
        end
    end
end
%日前市场发电机每一时间段出力上下限确定
P_DA_GTmax=zeros(size(AllGens,1),T);
P_DA_GTmin=zeros(size(AllGens,1),T);
for i=1:size(AllThermalGenerators,1)+size(AllWindFarmGenerators,1)+size(AllPhotovoltaicGenerators,1)
    for j=1:T
        if i<=size(AllThermalGenerators,1)
            P_DA_GTmax(i,j)=AllThermalGenerators{i}.capacity;
            P_DA_GTmin(i,j)=AllThermalGenerators{i}.P_min;
        elseif i>size(AllThermalGenerators,1)&&i<=size(AllThermalGenerators,1)+size(AllWindFarmGenerators,1)
            P_DA_GTmax(i,j)=AllWindFarmGenerators{i-size(AllThermalGenerators,1)}.typicaloutput(j);
            P_DA_GTmin(i,j)=0;
        else
            P_DA_GTmax(i,j)=AllPhotovoltaicGenerators{i-size(AllThermalGenerators,1)-size(AllWindFarmGenerators,1)}.typicaloutput(j);
            P_DA_GTmin(i,j)=0;
        end
    end
end
%实时市场发电机每报价段出力上下限确定
P_RT_Gmax=zeros(size(AllGens,1),Seg_offer,W,T);
P_RT_Gmin=zeros(size(AllGens,1),Seg_offer,W,T);
for i=1:size(AllThermalGenerators,1)+size(AllWindFarmGenerators,1)+size(AllPhotovoltaicGenerators,1)
    for j=1:Seg_offer
        for k=1:W
            for m=1:T
                if i<=size(AllThermalGenerators,1)
                    P_RT_Gmax(i,j,k,m)=AllThermalGenerators{i}.capacity;
                    P_RT_Gmin(i,j,k,m)=0;
                elseif i>size(AllThermalGenerators,1)&&i<=size(AllThermalGenerators,1)+size(AllWindFarmGenerators,1)
                    P_RT_Gmax(i,j,k,m)=AllWindFarmGenerators{i-size(AllThermalGenerators,1)}.typicaloutput(m)*AllScenarios{k}.WF_change;
                    P_RT_Gmin(i,j,k,m)=0;
                else
                    P_RT_Gmax(i,j,k,m)=AllPhotovoltaicGenerators{i-size(AllThermalGenerators,1)-size(AllWindFarmGenerators,1)}.typicaloutput(m)*AllScenarios{k}.PV_change;
                    P_RT_Gmin(i,j,k,m)=0;
                end
            end
        end
    end
end
%实时市场发电机每时间段出力上下限确定
P_RT_GTmax=zeros(size(AllGens,1),W,T);
P_RT_GTmin=zeros(size(AllGens,1),W,T);
for i=1:size(AllThermalGenerators,1)+size(AllWindFarmGenerators,1)+size(AllPhotovoltaicGenerators,1)
    for j=1:W
        for k=1:T
            if i<=size(AllThermalGenerators,1)
                P_RT_GTmax(i,j,k)=AllThermalGenerators{i}.capacity;
                P_RT_GTmin(i,j,k)=AllThermalGenerators{i}.P_min;
            elseif i>size(AllThermalGenerators,1)&&i<=size(AllThermalGenerators,1)+size(AllWindFarmGenerators,1)
                P_RT_GTmax(i,j,k)=AllWindFarmGenerators{i-size(AllThermalGenerators,1)}.typicaloutput(k)*AllScenarios{j}.WF_change;
                P_RT_GTmin(i,j,k)=0;
            else
                P_RT_GTmax(i,j,k)=AllPhotovoltaicGenerators{i-size(AllThermalGenerators,1)-size(AllWindFarmGenerators,1)}.typicaloutput(k)*AllScenarios{j}.PV_change;
                P_RT_GTmin(i,j,k)=0;
            end
        end
    end
end
%用电用户每一段报价上下限确定
P_Dmax=zeros(size(AllLoads,1),Seg_bid,T);
P_Dmin=zeros(size(AllLoads,1),Seg_bid,T);
for i=1:size(AllLoads,1)
    for j=1:Seg_bid
        for k=1:T
            P_Dmax(i,j,k)=AllLoads{i}.powerconsume_DA(k)/Seg_bid;
            P_Dmin(i,j,k)=0;
        end
    end
end
%用户的报价数据确定
Bid_Price=zeros(size(AllLoads,1), Seg_bid);
for i=1:size(AllLoads,1)
    for j=1:Seg_bid
        Bid_Price(i,j)=BidData{i+1,j+1};
    end
end
%TGC市场参数确定
Q_D=0;
for i=1:size(AllLoads,1)
    for j=1:Seg_bid
        for k=1:T
            Q_D=Q_D+P_Dmax(i,j,k);
        end
    end
end
Q_D=Q_D*T_Gap;
alpha_TGC=Penalty_Price;
beta_TGC=(1-Theta_FT)*Penalty_Price/(RPS_ratio*Q_D);
%实时市场收益线性化参数确定
P_max=zeros(size(RE_Gens,1),T); P_min=zeros(size(RE_Gens,1),T);
P_ref_gap=zeros(size(RE_Gens,1),T);
for i=1:size(RE_Gens,1)
    for j=1:T
        P_max(i,j)=P_DA_GTmax(RE_Gens(i),j);
        P_min(i,j)=P_DA_GTmin(RE_Gens(i),j);
        P_ref_gap(i,j)=(P_max(i,j)-P_min(i,j))/(Num_Approximate-1);
    end
end
disp('原始数据处理完毕');
%%
%求出一些基本数据
Num_Gen=size(AllGens,1);                                   %求出总发电机数目用于建模
Num_CoGen=size(RE_Gens,1);                                 %战略性可再生能源发电公司拥有的机组数目
Num_Load=size(AllLoads,1);                                 %求出总负荷数方便建模
Num_Node=size(AllNodes,1);                                 %求出总线路数用于建模
%创建变量
Offer_Price=sdpvar(Num_CoGen,Seg_offer,'full');              %创建发电机报价变量
LMP_DA=sdpvar(Num_Node,T,'full');                          %日前市场中的节点电价
P_DA_G=sdpvar(Num_Gen,Seg_offer,T,'full');                 %每一台发电机在日前市场不同时段以及报价块的出力
P_DA_D=sdpvar(Num_Load,Seg_bid,T,'full');                  %每一个用户在日前市场不同时段以及报价块的用电量
Delta_DA=sdpvar(Num_Node,T,'full');                        %日前市场每一个节点的电压角
LMP_RT=sdpvar(Num_Node,W,T,'full');                        %实时市场中的节点电价
P_RT_G=sdpvar(Num_Gen,Seg_offer,W,T,'full');               %每一台发电机在实时市场不同时段不同场景以及报价块的出力
P_RT_D=sdpvar(Num_Load,Seg_bid,W,T,'full');                %每一个用户在实时市场不同时段不同场景以及报价块的用电量
Delta_RT=sdpvar(Num_Node,W,T,'full');                      %实时市场每一个节点的电压角
TGC_Price_PT=sdpvar(W,1,'full');                           %当前市场中的TGC价格
TGC_Price_FT=sdpvar(1,1,'full');                           %未来市场中的TGC价格
Q_PT=sdpvar(W,1,'full');                                   %在当前市场中出售的TGC数量
Q_PT_Other=sdpvar(W,1,'full');                             %其他可再生能源发电公司在当前市场出售的TGC数量
Q_FT=sdpvar(W,1,'full');                                   %在未来市场中出售的TGC数量
Q_Gmax=sdpvar(W,1,'full');                                 %战略性可再生能源发电公司当前市场中出售的最大TGC数量
Q_Gmax_Other=sdpvar(W,1,'full');                           %其他可再生能源发电公司当前市场中出售的最大TGC数量
%创建线性化过程中产生的变量
Z_DA1=sdpvar(T,1,'full');                                  %日前市场中的线性化项1
Z_DA2=sdpvar(T,1,'full');                                  %日前市场中的线性化项2
Z_RT11=sdpvar(W,T,'full');                                 %实时市场中的线性化项1-1
Z_RT12=sdpvar(W,T,'full');                                 %实时市场中的线性化项1-2
Z_RT2=sdpvar(W,T,'full');                                  %实时市场中的线性化项2
Z_TGC=sdpvar(W,1,'full');                                  %TGC市场中的线性化项
P_ref=sdpvar(Num_CoGen,Num_Approximate,'full');            %实时市场线性化中的近似基准值
RT_Ancillary1=sdpvar(Num_CoGen,W,T,Num_Approximate,'full');%RT线性化辅助变量1
RT_Ancillary2=binvar(Num_CoGen,T,Num_Approximate,'full');  %RT线性化辅助变量2
Q_ref=sdpvar(Num_Approximate,1,'full');                    %TGC市场线性化中的近似基准值
TGC_Ancillary1=sdpvar(W,Num_Approximate,'full');           %TGC线性化辅助变量1
TGC_Ancillary2=binvar(W,Num_Approximate,'full');           %TGC线性化辅助变量2
%创建KKT条件中的拉格朗日乘子变量
u_DA_Gmin=sdpvar(Num_Gen,Seg_offer,T,'full');              %日前市场发电机每一时段每一报价段出力约束乘子
u_DA_Gmax=sdpvar(Num_Gen,Seg_offer,T,'full');
u_DA_GTmin=sdpvar(Num_Gen,T,'full');                       %日前市场发电机每一时段出力约束乘子
u_DA_GTmax=sdpvar(Num_Gen,T,'full');
u_DA_Dmin=sdpvar(Num_Load,Seg_bid,T,'full');               %日前市场用户每一时段每一用电段用电量约束乘子
u_DA_Dmax=sdpvar(Num_Load,Seg_bid,T,'full');
v_DA_Lmax=sdpvar(Num_Node,Num_Node,T,'full');              %日前市场线路流过潮流约束乘子
e_DA_min=sdpvar(Num_Node,T,'full');                        %日前市场电压角约束乘子
e_DA_max=sdpvar(Num_Node,T,'full');
e_DA1=sdpvar(T,1,'full');                                  %日前市场参考节点电压角约束乘子
u_RT_Gmin=sdpvar(Num_Gen,Seg_offer,W,T,'full');            %实时市场发电机每一时段每一报价段每一场景出力约束乘子
u_RT_Gmax=sdpvar(Num_Gen,Seg_offer,W,T,'full');
u_RT_GTmin=sdpvar(Num_Gen,W,T,'full');                     %实时市场发电机每一时段每一场景出力约束乘子
u_RT_GTmax=sdpvar(Num_Gen,W,T,'full');
u_RT_Dmin=sdpvar(Num_Load,Seg_bid,W,T,'full');             %实时市场用户每一时段每一用电段每一场景用电量约束乘子
u_RT_Dmax=sdpvar(Num_Load,Seg_bid,W,T,'full');
v_RT_Lmax=sdpvar(Num_Node,Num_Node,W,T,'full');            %实时市场线路流过潮流约束乘子
e_RT_min=sdpvar(Num_Node,W,T,'full');                      %实时市场电压角约束乘子
e_RT_max=sdpvar(Num_Node,W,T,'full');
e_RT1=sdpvar(W,T,'full'); 
%创建KKT互补松弛约束线性化中的变量
c_DA_Gmin=binvar(Num_Gen,Seg_offer,T,'full');              %日前市场发电机每一时段每一报价段出力约束线性化0-1变量
c_DA_Gmax=binvar(Num_Gen,Seg_offer,T,'full');
c_DA_GTmin=binvar(Num_Gen,T,'full');                       %日前市场发电机每一时段出力约束线性化0-1变量
c_DA_GTmax=binvar(Num_Gen,T,'full');
c_DA_Dmin=binvar(Num_Load,Seg_bid,T,'full');               %日前市场用户每一时段每一用电段用电量约束线性化0-1变量
c_DA_Dmax=binvar(Num_Load,Seg_bid,T,'full');
c_DA_Lmax=binvar(Num_Node,Num_Node,T,'full');              %日前市场线路流过潮流约束线性化0-1变量
c_DA_Jmin=binvar(Num_Node,T,'full');                       %日前市场电压角约束线性化0-1变量
c_DA_Jmax=binvar(Num_Node,T,'full');
c_RT_Gmin=binvar(Num_Gen,Seg_offer,W,T,'full');            %实时市场发电机每一时段每一报价段每一场景出力约束线性化0-1变量
c_RT_Gmax=binvar(Num_Gen,Seg_offer,W,T,'full');
c_RT_GTmin=binvar(Num_Gen,W,T,'full');                     %实时市场发电机每一时段每一场景出力约束线性化0-1变量
c_RT_GTmax=binvar(Num_Gen,W,T,'full');
c_RT_Dmin=binvar(Num_Load,Seg_bid,W,T,'full');             %实时市场用户每一时段每一用电段每一场景用电量约束线性化0-1变量
c_RT_Dmax=binvar(Num_Load,Seg_bid,W,T,'full');
c_RT_Lmax=binvar(Num_Node,Num_Node,W,T,'full');            %实时市场线路流过潮流约束线性化0-1变量
c_RT_Jmin=binvar(Num_Node,W,T,'full');                     %实时市场电压角约束线性化0-1变量
c_RT_Jmax=binvar(Num_Node,W,T,'full');
%添加松弛
SC=sdpvar(1,1,'full');
disp('变量创建完毕');
%%
%给初始值
load('Initial_RT.mat');
assign(c_RT_Dmax,c_RT_Dmax_Ini);
assign(c_RT_Dmin,c_RT_Dmin_Ini);
assign(c_RT_Gmax,c_RT_Gmax_Ini);
assign(c_RT_Gmin,c_RT_Gmin_Ini);
assign(c_RT_GTmax,c_RT_GTmax_Ini);
assign(c_RT_GTmin,c_RT_GTmin_Ini);
assign(c_RT_Jmax,c_RT_Jmax_Ini);
assign(c_RT_Jmin,c_RT_Jmin_Ini);
assign(c_RT_Lmax,c_RT_Lmax_Ini);
assign(Delta_RT,Delta_RT_Ini);
assign(e_RT1,e_RT1_Ini);
assign(e_RT_max,e_RT_max_Ini);
assign(e_RT_min,e_RT_min_Ini);
assign(LMP_RT,LMP_RT_Ini);
assign(P_RT_D,P_RT_D_Ini);
assign(P_RT_G,P_RT_G_Ini);
assign(u_RT_Dmax,u_RT_Dmax_Ini);
assign(u_RT_Dmin,u_RT_Dmin_Ini);
assign(u_RT_Gmax,u_RT_Gmax_Ini);
assign(u_RT_Gmin,u_RT_Gmin_Ini);
assign(u_RT_GTmax,u_RT_GTmax_Ini);
assign(u_RT_GTmin,u_RT_GTmin_Ini);
assign(v_RT_Lmax,v_RT_Lmax_Ini);
%%
%添加约束
Const_Upper=[];
%顶层问题约束     
for i=1:Num_CoGen                                         %报价上下限约束            
    for j=1:Seg_offer
        Const_Upper=[Const_Upper,Offer_Price(i,j)<=Offer_UpperCup];
        Const_Upper=[Const_Upper,Offer_Price(i,j)>=Offer_LowerCup];
    end
end

for i=1:Num_CoGen                                         %报价块递增约束            
    for j=2:Seg_offer
        Const_Upper=[Const_Upper,Offer_Price(i,j-1)<=Offer_Price(i,j)];
    end
end

for i=1:W                                                %TGC出售数量约束(下限暂时设为0)    
    Const_Upper=[Const_Upper,Q_PT(i)<=Q_Gmax(i)];
    Const_Upper=[Const_Upper,Q_PT(i)>=0];
end

for i=1:W                                                %战略性可再生能源发电公司TGC出售数量上限约束
    New_Const=0;
    for j=1:size(RE_Gens,1)
        for k=1:Seg_offer
            for m=1:T
                New_Const=New_Const+P_RT_G(RE_Gens(j),k,i,m);
            end
        end
    end
    Const_Upper=[Const_Upper,Q_Gmax(i)==New_Const*T_Gap];
end

for i=1:W                                                %其他可再生能源发电公司TGC出售数量上限约束
    New_Const=0;
    for j=1:size(RE_Gens_Other,1)
        for k=1:Seg_offer
            for m=1:T
                New_Const=New_Const+P_RT_G(RE_Gens_Other(j),k,i,m);
            end
        end
    end
    Const_Upper=[Const_Upper,Q_Gmax_Other(i)==New_Const*T_Gap];
end

for i=1:W                                                %Q_FT计算公式约束
    Const_Upper=[Const_Upper,Q_FT(i)==Q_Gmax(i)-Q_PT(i)];
end

for i=1:W                                                %Q_PT_Other计算公式约束
    Const_Upper=[Const_Upper,Q_PT_Other(i)==Q_Gmax_Other(i)];
end

disp('顶层问题约束添加完毕');
%%
%底层问题一：日前市场出清约束
Const_DA=[];

for i=1:Num_Gen                                          %KKT条件转换后约束1
    for j=1:Seg_offer
        for k=1:T
            if ismember(i,RE_Gens)
                Const_DA=[Const_DA,Offer_Price(RE_Gens==i,j)-LMP_DA(AllGens{i}.location,k)+u_DA_Gmax(i,j,k)-u_DA_Gmin(i,j,k)+u_DA_GTmax(i,k)-u_DA_GTmin(i,k)==0];
            else
                Const_DA=[Const_DA,Marign_Cost(i,j)-LMP_DA(AllGens{i}.location,k)+u_DA_Gmax(i,j,k)-u_DA_Gmin(i,j,k)+u_DA_GTmax(i,k)-u_DA_GTmin(i,k)==0];
            end
        end
    end
end

for i=1:Num_Load                                         %KKT条件转换后约束2
    for j=1:Seg_bid
        for k=1:T
            Const_DA=[Const_DA,LMP_DA(AllLoads{i}.location,k)-Bid_Price(i,j)+u_DA_Dmax(i,j,k)-u_DA_Dmin(i,j,k)==0];
        end
    end
end

for i=1:Num_Node                                         %KKT条件转换后约束3
    for j=1:T
        New_Const=0;
        for k=1:size(AllNodes{i}.connectnodes,1)
            New_Const=New_Const+B(i,AllNodes{i}.connectnodes(k))*(LMP_DA(i,j)-LMP_DA(AllNodes{i}.connectnodes(k),j))+B(i,AllNodes{i}.connectnodes(k))*(v_DA_Lmax(i,AllNodes{i}.connectnodes(k),j)-v_DA_Lmax(AllNodes{i}.connectnodes(k),i,j));
        end
        if i~=1
            Const_DA=[Const_DA,New_Const+e_DA_max(i,j)-e_DA_min(i,j)==0];
        else
            Const_DA=[Const_DA,New_Const+e_DA_max(i,j)-e_DA_min(i,j)+e_DA1(j)==0];
        end      
    end
end

for i=1:Num_Gen                                         %日前市场发电机每一时段每一报价段出力互补松弛约束
    for j=1:Seg_offer
        for k=1:T
            Const_DA=[Const_DA,0<=(P_DA_G(i,j,k)-P_DA_Gmin(i,j,k))<=M*c_DA_Gmin(i,j,k)];
            Const_DA=[Const_DA,0<=u_DA_Gmin(i,j,k)<=M*(1-c_DA_Gmin(i,j,k))];
            Const_DA=[Const_DA,0<=(P_DA_Gmax(i,j,k)-P_DA_G(i,j,k))<=M*c_DA_Gmax(i,j,k)];
            Const_DA=[Const_DA,0<=u_DA_Gmax(i,j,k)<=M*(1-c_DA_Gmax(i,j,k))];
        end
    end
end

for i=1:Num_Gen                                         %日前市场发电机每一时段出力互补松弛约束
    for j=1:T
        Const_DA=[Const_DA,0<=(sum(P_DA_G(i,:,j))-P_DA_GTmin(i,j))<=M*c_DA_GTmin(i,j)];
        Const_DA=[Const_DA,0<=u_DA_GTmin(i,j)<=M*(1-c_DA_GTmin(i,j))];
        Const_DA=[Const_DA,0<=(P_DA_GTmax(i,j)-sum(P_DA_G(i,:,j)))<=M*c_DA_GTmax(i,j)];
        Const_DA=[Const_DA,0<=u_DA_GTmax(i,j)<=M*(1-c_DA_GTmax(i,j))];
    end
end

for i=1:Num_Load                                        %日前市场用户每一时段每一用电段用电量互补松弛约束
    for j=1:Seg_bid
        for k=1:T
            Const_DA=[Const_DA,0<=(P_DA_D(i,j,k)-P_Dmin(i,j,k))<=M*c_DA_Dmin(i,j,k)];
            Const_DA=[Const_DA,0<=u_DA_Dmin(i,j,k)<=M*(1-c_DA_Dmin(i,j,k))];
            Const_DA=[Const_DA,0<=(P_Dmax(i,j,k)-P_DA_D(i,j,k))<=M*c_DA_Dmax(i,j,k)];
            Const_DA=[Const_DA,0<=u_DA_Dmax(i,j,k)<=M*(1-c_DA_Dmax(i,j,k))];
        end
    end
end

for i=1:Num_Node                                        %日前市场线路流过潮流互补松弛约束
    for j=1:size(AllNodes{i}.connectnodes,1)
        for k=1:T
            Const_DA=[Const_DA,0<=P_Lmax(i,AllNodes{i}.connectnodes(j))-B(i,AllNodes{i}.connectnodes(j))*(Delta_DA(i,k)-Delta_DA(AllNodes{i}.connectnodes(j),k))<=M*c_DA_Lmax(i,AllNodes{i}.connectnodes(j),k)];
            Const_DA=[Const_DA,0<=v_DA_Lmax(i,AllNodes{i}.connectnodes(j),k)<=M*(1-c_DA_Lmax(i,AllNodes{i}.connectnodes(j),k))];
        end
    end
end

for i=1:Num_Node                                        %日前市场电压角互补松弛约束
    for j=1:T
        Const_DA=[Const_DA,0<=Delta_DA(i,j)+pi<=M*c_DA_Jmin(i,j)];
        Const_DA=[Const_DA,0<=e_DA_min(i,j)<=M*(1-c_DA_Jmin(i,j))];
        Const_DA=[Const_DA,0<=pi-Delta_DA(i,j)<=M*c_DA_Jmax(i,j)];
        Const_DA=[Const_DA,0<=e_DA_max(i,j)<=M*(1-c_DA_Jmax(i,j))];
    end
end

for i=1:Num_Node                                        %日前市场功率平衡约束
    for j=1:T
        New_Const=0;
        for k=1:Num_Gen
            if AllGens{k}.location==i
                New_Const=New_Const+sum(P_DA_G(k,:,j));
            end
        end
        for k=1:Num_Load
            if AllLoads{k}.location==i
                New_Const=New_Const-sum(P_DA_D(k,:,j));
            end
        end
        for k=1:size(AllNodes{i}.connectnodes,1)
            New_Const=New_Const-B(i,AllNodes{i}.connectnodes(k))*(Delta_DA(i,j)-Delta_DA(AllNodes{i}.connectnodes(k),j));
        end
        Const_DA=[Const_DA,New_Const==0];
    end
end

for i=1:T                                              %日前市场参考点电压角约束
    Const_DA=[Const_DA,Delta_DA(1,i)==0];
end
disp('底层问题一问题约束添加完毕');
%%
%底层问题二：实时市场出清约束
Const_RT=[];
for m=1:W
    for i=1:Num_Gen                                          %KKT条件转换后约束1
        for j=1:Seg_offer
            for k=1:T
                if ismember(i,RE_Gens)
                    Const_RT=[Const_RT,Offer_Price(RE_Gens==i,j)-LMP_RT(AllGens{i}.location,m,k)+u_RT_Gmax(i,j,m,k)-u_RT_Gmin(i,j,m,k)+u_RT_GTmax(i,m,k)-u_RT_GTmin(i,m,k)==0];
                else
                    Const_RT=[Const_RT,Marign_Cost(i,j)-LMP_RT(AllGens{i}.location,m,k)+u_RT_Gmax(i,j,m,k)-u_RT_Gmin(i,j,m,k)+u_RT_GTmax(i,m,k)-u_RT_GTmin(i,m,k)==0];
                end
            end
        end
    end

    for i=1:Num_Load                                         %KKT条件转换后约束2
        for j=1:Seg_bid
            for k=1:T
                Const_RT=[Const_RT,LMP_RT(AllLoads{i}.location,m,k)-Bid_Price(i,j)+u_RT_Dmax(i,j,m,k)-u_RT_Dmin(i,j,m,k)==0];
            end
        end
    end

    for i=1:Num_Node                                         %KKT条件转换后约束3
        for j=1:T
            New_Const=0;
            for k=1:size(AllNodes{i}.connectnodes,1)
                New_Const=New_Const+B(i,AllNodes{i}.connectnodes(k))*(LMP_RT(i,m,j)-LMP_RT(AllNodes{i}.connectnodes(k),m,j))+B(i,AllNodes{i}.connectnodes(k))*(v_RT_Lmax(i,AllNodes{i}.connectnodes(k),m,j)-v_RT_Lmax(AllNodes{i}.connectnodes(k),i,m,j));
            end
            if i~=1
                Const_RT=[Const_RT,New_Const+e_RT_max(i,m,j)-e_RT_min(i,m,j)==0];
            else
                Const_RT=[Const_RT,New_Const+e_RT_max(i,m,j)-e_RT_min(i,m,j)+e_RT1(m,j)==0];
            end      
        end
    end

    for i=1:Num_Gen                                         %实时市场发电机每一时段每一报价段出力互补松弛约束
        for j=1:Seg_offer
            for k=1:T
                Const_RT=[Const_RT,0<=(P_RT_G(i,j,m,k)-P_RT_Gmin(i,j,m,k))<=M*c_RT_Gmin(i,j,m,k)];
                Const_RT=[Const_RT,0<=u_RT_Gmin(i,j,m,k)<=M*(1-c_RT_Gmin(i,j,m,k))];
                Const_RT=[Const_RT,0<=(P_RT_Gmax(i,j,m,k)-P_RT_G(i,j,m,k))<=M*c_RT_Gmax(i,j,m,k)];
                Const_RT=[Const_RT,0<=u_RT_Gmax(i,j,m,k)<=M*(1-c_RT_Gmax(i,j,m,k))];
            end
        end
    end

    for i=1:Num_Gen                                         %实时市场发电机每一时段出力互补松弛约束
        for j=1:T
            Const_RT=[Const_RT,0<=(sum(P_RT_G(i,:,m,j))-P_RT_GTmin(i,m,j))<=M*c_RT_GTmin(i,m,j)];
            Const_RT=[Const_RT,0<=u_RT_GTmin(i,m,j)<=M*(1-c_RT_GTmin(i,m,j))];
            Const_RT=[Const_RT,0<=(P_RT_GTmax(i,m,j)-sum(P_RT_G(i,:,m,j)))<=M*c_RT_GTmax(i,m,j)];
            Const_RT=[Const_RT,0<=u_RT_GTmax(i,m,j)<=M*(1-c_RT_GTmax(i,m,j))];
        end
    end

    for i=1:Num_Load                                        %实时市场用户每一时段每一用电段用电量互补松弛约束
        for j=1:Seg_bid
            for k=1:T
                Const_RT=[Const_RT,0<=(P_RT_D(i,j,m,k)-P_Dmin(i,j,k))<=M*c_RT_Dmin(i,j,m,k)];
                Const_RT=[Const_RT,0<=u_RT_Dmin(i,j,m,k)<=M*(1-c_RT_Dmin(i,j,m,k))];
                Const_RT=[Const_RT,0<=(P_Dmax(i,j,k)-P_RT_D(i,j,m,k))<=M*c_RT_Dmax(i,j,m,k)];
                Const_RT=[Const_RT,0<=u_RT_Dmax(i,j,m,k)<=M*(1-c_RT_Dmax(i,j,m,k))];
            end
        end
    end

    for i=1:Num_Node                                        %实时市场线路流过潮流互补松弛约束
        for j=1:size(AllNodes{i}.connectnodes,1)
            for k=1:T
                Const_RT=[Const_RT,0<=P_Lmax(i,AllNodes{i}.connectnodes(j))-B(i,AllNodes{i}.connectnodes(j))*(Delta_RT(i,m,k)-Delta_RT(AllNodes{i}.connectnodes(j),m,k))<=M*c_RT_Lmax(i,AllNodes{i}.connectnodes(j),m,k)];
                Const_RT=[Const_RT,0<=v_RT_Lmax(i,AllNodes{i}.connectnodes(j),m,k)<=M*(1-c_RT_Lmax(i,AllNodes{i}.connectnodes(j),m,k))];
            end
        end
    end

    for i=1:Num_Node                                        %实时市场电压角互补松弛约束
        for j=1:T
            Const_RT=[Const_RT,0<=Delta_RT(i,m,j)+pi<=M*c_RT_Jmin(i,m,j)];
            Const_RT=[Const_RT,0<=e_RT_min(i,m,j)<=M*(1-c_RT_Jmin(i,m,j))];
            Const_RT=[Const_RT,0<=pi-Delta_RT(i,m,j)<=M*c_RT_Jmax(i,m,j)];
            Const_RT=[Const_RT,0<=e_RT_max(i,m,j)<=M*(1-c_RT_Jmax(i,m,j))];
        end
    end

    for i=1:Num_Node                                        %实时市场功率平衡约束
        for j=1:T
            New_Const=0;
            for k=1:Num_Gen
                if AllGens{k}.location==i
                    New_Const=New_Const+sum(P_RT_G(k,:,m,j));
                end
            end
            for k=1:Num_Load
                if AllLoads{k}.location==i
                    New_Const=New_Const-sum(P_RT_D(k,:,m,j));
                end
            end
            for k=1:size(AllNodes{i}.connectnodes,1)
                New_Const=New_Const-B(i,AllNodes{i}.connectnodes(k))*(Delta_RT(i,m,j)-Delta_RT(AllNodes{i}.connectnodes(k),m,j));
            end
            Const_RT=[Const_RT,New_Const==0];
        end
    end

    for i=1:T                                              %实时市场参考点电压角约束
        Const_RT=[Const_RT,Delta_RT(1,m,i)==0];
    end
end
disp('底层问题二问题约束添加完毕');
%%
Const_Liner_DA=[];

%日前市场收益线性化约束
for i=1:T                                            %日前市场线性化项Z_DA1计算公式
    New_Const=0;
    for j=1:Num_Load
        for k=1:Seg_bid
            New_Const=New_Const+Bid_Price(j,k)*P_DA_D(j,k,i);
        end
    end
    for j=1:Num_Gen
        if ~ismember(j,RE_Gens)
            for k=1:Seg_offer
                New_Const=New_Const-Marign_Cost(j,k)*P_DA_G(j,k,i);
            end
        end
    end
    for j=1:Num_Gen
        for k=1:Seg_offer
            New_Const=New_Const+u_DA_Gmin(j,k,i)*P_DA_Gmin(j,k,i)-u_DA_Gmax(j,k,i)*P_DA_Gmax(j,k,i);
        end
    end
    for j=1:Num_Gen
        New_Const=New_Const+u_DA_GTmin(j,i)*P_DA_GTmin(j,i)-u_DA_GTmax(j,i)*P_DA_GTmax(j,i);
    end
    for j=1:Num_Load
        for k=1:Seg_bid
            New_Const=New_Const+u_DA_Dmin(j,k,i)*P_Dmin(j,k,i)-u_DA_Dmax(j,k,i)*P_Dmax(j,k,i);
        end
    end
    for j=1:Num_Node
        for k=1:size(AllNodes{j}.connectnodes,1)
            New_Const=New_Const-v_DA_Lmax(j,AllNodes{j}.connectnodes(k),i)*P_Lmax(j,AllNodes{j}.connectnodes(k));
        end
    end
    for j=1:Num_Node
        New_Const=New_Const-e_DA_max(j,i)*pi-e_DA_min(j,i)*pi;
    end
    Const_Liner_DA=[Const_Liner_DA,Z_DA1(i)==New_Const];
end

for i=1:T                                            %日前市场线性化项Z_DA2计算公式
    New_Const=0;
    for j=1:size(RE_Gens,1)
        for k=1:Seg_offer
            New_Const=New_Const+u_DA_Gmax(RE_Gens(j),k,i)*P_DA_Gmax(RE_Gens(j),k,i)-u_DA_Gmin(RE_Gens(j),k,i)*P_DA_Gmin(RE_Gens(j),k,i)-Marign_Cost(RE_Gens(j),k)*P_DA_G(RE_Gens(j),k,i);
        end
    end
    for j=1:size(RE_Gens,1)
        New_Const=New_Const+u_DA_GTmax(RE_Gens(j),i)*P_DA_GTmax(RE_Gens(j),i)-u_DA_GTmin(RE_Gens(j),i)*P_DA_GTmin(RE_Gens(j),i);
    end
    Const_Liner_DA=[Const_Liner_DA,Z_DA2(i)==New_Const];
end
disp('日前市场线性化项约束添加完毕');
%%
Const_Liner_RT=[];

%实时市场收益线性化约束
for i=1:W                                            %实时市场线性化项Z_RT11计算公式
    for j=1:T
        New_Const=0;
        for k=1:Num_Load
            for m=1:Seg_bid
                New_Const=New_Const+Bid_Price(k,m)*P_RT_D(k,m,i,j);
            end
        end
        for k=1:Num_Gen
            for m=1:Seg_offer
                if ~ismember(k,RE_Gens)
                    New_Const=New_Const-Marign_Cost(k,m)*P_RT_G(k,m,i,j);
                end
            end
        end
        for k=1:Num_Gen
            for m=1:Seg_offer
                New_Const=New_Const+u_RT_Gmin(k,m,i,j)*P_RT_Gmin(k,m,i,j)-u_RT_Gmax(k,m,i,j)*P_RT_Gmax(k,m,i,j);
            end
        end
        for k=1:Num_Gen
            New_Const=New_Const+u_RT_GTmin(k,i,j)*P_RT_GTmin(k,i,j)-u_RT_GTmax(k,i,j)*P_RT_GTmax(k,i,j);
        end
        for k=1:Num_Load
            for m=1:Seg_bid
                New_Const=New_Const+u_RT_Dmin(k,m,i,j)*P_Dmin(k,m,j)-u_RT_Dmax(k,m,i,j)*P_Dmax(k,m,j);
            end
        end
        for k=1:Num_Node
            for m=1:size(AllNodes{k}.connectnodes,1)
                New_Const=New_Const-v_RT_Lmax(k,AllNodes{k}.connectnodes(m),i,j)*P_Lmax(k,AllNodes{k}.connectnodes(m));
            end
        end
        for k=1:Num_Node
            New_Const=New_Const-e_RT_min(k,i,j)*pi-e_RT_max(k,i,j)*pi;
        end
        Const_Liner_RT=[Const_Liner_RT,Z_RT11(i,j)==New_Const];
    end
end

for i=1:W                                            %实时市场线性化项Z_RT12计算公式        
    for j=1:T
        New_Const=0;
        for k=1:size(RE_Gens,1)
            for m=1:Num_Approximate
                New_Const=New_Const-(P_min(k,j)+(m-1)*P_ref_gap(k,j))*RT_Ancillary1(k,i,j,m);
            end
        end
        Const_Liner_RT=[Const_Liner_RT,Z_RT12(i,j)==New_Const];
    end
end
for i=1:size(RE_Gens,1)
    for j=1:T
        New_Const=0;
        for k=1:Num_Approximate
            New_Const=New_Const+(P_min(i,j)+(k-1)*P_ref_gap(i,j))*RT_Ancillary2(i,j,k);
        end
        Const_Liner_RT=[Const_Liner_RT,sum(P_DA_G(RE_Gens(i),:,j))-P_ref_gap/2<=New_Const];
        Const_Liner_RT=[Const_Liner_RT,sum(P_DA_G(RE_Gens(i),:,j))+P_ref_gap/2>=New_Const];
    end
end
for i=1:size(RE_Gens,1)
    for j=1:W
        for k=1:T
            for m=1:Num_Approximate
                Const_Liner_RT=[Const_Liner_RT,0<=LMP_RT(AllGens{RE_Gens(i)}.location,j,k)-RT_Ancillary1(i,j,k,m)<=M*(1-RT_Ancillary2(i,k,m))];
                Const_Liner_RT=[Const_Liner_RT,0<=RT_Ancillary1(i,j,k,m)<=M*RT_Ancillary2(i,k,m)];
            end
        end
    end
end
for i=1:size(RE_Gens,1)
    for j=1:T
        Const_Liner_RT=[Const_Liner_RT,sum(RT_Ancillary2(i,j,:))==1];
    end
end

for i=1:W                                            %实时市场线性化项Z_RT2计算公式        
    for j=1:T
        New_Const=0;
        for k=1:size(RE_Gens,1)
            for m=1:Seg_offer
                New_Const=New_Const-Marign_Cost(RE_Gens(k),m)*(P_RT_G(RE_Gens(k),m,i,j)-P_DA_G(RE_Gens(k),m,j));
            end
        end
        Const_Liner_RT=[Const_Liner_RT,Z_RT2(i,j)==New_Const];
    end
end
disp('实时市场线性化项约束添加完毕');
%%
Const_Liner_TGC=[];

%TGC市场收益线性化约束
for i=1:W                                           %TGC市场Z_TGC计算公式
    New_Const=0;
    for j=1:Num_Approximate
        New_Const=New_Const+(Q_min+(j-1)*Q_ref_gap)*TGC_Ancillary1(i,j);
    end
    Const_Liner_TGC=[Const_Liner_TGC,Z_TGC(i)==New_Const];
end
for i=1:W
    New_Const=0;
    for j=1:Num_Approximate
        New_Const=New_Const+(Q_min+(j-1)*Q_ref_gap)*TGC_Ancillary2(i,j);
    end
    Const_Liner_TGC=[Const_Liner_TGC,Q_PT(i)-Q_ref_gap/2<=New_Const];
    Const_Liner_TGC=[Const_Liner_TGC,Q_PT(i)+Q_ref_gap/2>=New_Const];
end
for i=1:W
    for j=1:Num_Approximate
        Const_Liner_TGC=[Const_Liner_TGC,0<=Q_PT_Other(i)-TGC_Ancillary1(i,j)<=M*(1-TGC_Ancillary2(i,j))];
        Const_Liner_TGC=[Const_Liner_TGC,0<=TGC_Ancillary1(i,j)<=M*TGC_Ancillary2(i,j)];
    end
end
for i=1:W
    Const_Liner_TGC=[Const_Liner_TGC,sum(TGC_Ancillary2(i,:))==1];
end
disp('TGC市场线性化项约束添加完毕');
%%
%计算目标函数
R_DA=(sum(Z_DA1)+sum(Z_DA2))*T_Gap;                       %日前市场收益计算 

New_Const=0;                                              %实时市场收益计算      
for i=1:W
    New_Const=New_Const+AllScenarios{i}.probability*(sum(Z_RT11(i,:))+sum(Z_RT12(i,:))+sum(Z_RT2(i,:)))*T_Gap;   
end
R_RT=New_Const;

New_Const=0;                                              %TGC市场收益计算约束       
for i=1:W
    New_Const=New_Const+AllScenarios{i}.probability*(alpha_TGC*Q_PT(i)-beta_TGC*Z_TGC(i)-beta_TGC*Q_PT(i)*Q_PT(i)+Gama_FT*Penalty_Price*(Q_Gmax(i)-Q_PT(i)));   
end
R_TGC=New_Const;

Obj=-(R_DA+R_RT+R_TGC);
%开始求解
Ops=sdpsettings('verbose',1,'solver','gurobi','usex0',1);    %设置求解参数
%Const=[Const_Upper,Const_DA,Const_RT,Const_Liner_DA,Const_Liner_RT,Const_Liner_TGC];
disp('开始求解');

Result=solvesdp([Const_Upper,Const_RT],0,Ops);
if Result.problem==0 
    fprintf('求解成功，输出信息为：%s\n',Result.info);
else
    fprintf('求解过程中出错，错误信息为：%s\n',Result.info);
end
%%
%日前市场约束正确性检验
Num_Const_DA=length(Const_DA);
Check_DA=cell(Num_Const_DA,2);
now=0; Iswrong_DA=0;
for i=1:Num_Gen                                          %KKT条件转换后约束1
    for j=1:Seg_offer
        for k=1:T
            now=now+1;
            if ismember(i,RE_Gens)
                Check_DA{now,1}=(abs(value(Offer_Price(RE_Gens==i,j)-LMP_DA(AllGens{i}.location,k)+u_DA_Gmax(i,j,k)-u_DA_Gmin(i,j,k)+u_DA_GTmax(i,k)-u_DA_GTmin(i,k)))<=gap);
                Check_DA{now,2}='KKT条件转换后约束1';
            else
                Check_DA{now,1}=(abs(value(Marign_Cost(i,j)-LMP_DA(AllGens{i}.location,k)+u_DA_Gmax(i,j,k)-u_DA_Gmin(i,j,k)+u_DA_GTmax(i,k)-u_DA_GTmin(i,k)))<=gap);
                Check_DA{now,2}='KKT条件转换后约束1';
            end
        end
    end
end

for i=1:Num_Load                                         %KKT条件转换后约束2
    for j=1:Seg_bid
        for k=1:T
            now=now+1;
            Check_DA{now,1}=(abs(value(LMP_DA(AllLoads{i}.location,k)-Bid_Price(i,j)+u_DA_Dmax(i,j,k)-u_DA_Dmin(i,j,k)))<=gap);
            Check_DA{now,2}='KKT条件转换后约束2';
        end
    end
end

for i=1:Num_Node                                         %KKT条件转换后约束3
    for j=1:T
        now=now+1;
        New_Const=0;
        for k=1:size(AllNodes{i}.connectnodes,1)
            New_Const=New_Const+B(i,AllNodes{i}.connectnodes(k))*(LMP_DA(i,j)-LMP_DA(AllNodes{i}.connectnodes(k),j))+B(i,AllNodes{i}.connectnodes(k))*(v_DA_Lmax(i,AllNodes{i}.connectnodes(k),j)-v_DA_Lmax(AllNodes{i}.connectnodes(k),i,j));
        end
        if i~=1
            Check_DA{now,1}=(abs(value(New_Const+e_DA_max(i,j)-e_DA_min(i,j)))<=gap);
            Check_DA{now,2}='KKT条件转换后约束3';
        else
            Check_DA{now,1}=(abs(value(New_Const+e_DA_max(i,j)-e_DA_min(i,j)+e_DA1(j)))<=gap);
            Check_DA{now,2}='KKT条件转换后约束3';
        end      
    end
end

for i=1:Num_Gen                                         %日前市场发电机每一时段每一报价段出力互补松弛约束
    for j=1:Seg_offer
        for k=1:T
            now=now+1;
            Check_DA{now,1}=value(-gap<=(P_DA_G(i,j,k)-P_DA_Gmin(i,j,k))<=M*c_DA_Gmin(i,j,k)+gap);
            Check_DA{now,2}='日前市场发电机出力下限约束';
            
            now=now+1;
            Check_DA{now,1}=value(-gap<=u_DA_Gmin(i,j,k)<=M*(1-c_DA_Gmin(i,j,k))+gap);
            Check_DA{now,2}='日前市场发电机出力下限拉格朗日乘子约束';
            
            now=now+1;
            Check_DA{now,1}=value(-gap<=(P_DA_Gmax(i,j,k)-P_DA_G(i,j,k))<=M*c_DA_Gmax(i,j,k)+gap);
            Check_DA{now,2}='日前市场发电机出力上限约束';
            
            now=now+1;
            Check_DA{now,1}=value(-gap<=u_DA_Gmax(i,j,k)<=M*(1-c_DA_Gmax(i,j,k))+gap);
            Check_DA{now,2}='日前市场发电机出力上限拉格朗日乘子约束';
        end
    end
end

for i=1:Num_Gen                                         %日前市场发电机每一时段出力互补松弛约束
    for j=1:T
        now=now+1;
        Check_DA{now,1}=value(-gap<=(sum(P_DA_G(i,:,j))-P_DA_GTmin(i,j))<=M*c_DA_GTmin(i,j)+gap);
        Check_DA{now,2}='日前市场发电机总出力下限约束';
        
        now=now+1;
        Check_DA{now,1}=value(-gap<=u_DA_GTmin(i,j)<=M*(1-c_DA_GTmin(i,j))+gap);
        Check_DA{now,2}='日前市场发电机总出力下限拉格朗日乘子约束';
        
        now=now+1;
        Check_DA{now,1}=value(-gap<=(P_DA_GTmax(i,j)-sum(P_DA_G(i,:,j)))<=M*c_DA_GTmax(i,j)+gap);
        Check_DA{now,2}='日前市场发电机总出力上限约束';
        
        now=now+1;
        Check_DA{now,1}=value(-gap<=u_DA_GTmax(i,j)<=M*(1-c_DA_GTmax(i,j))+gap);
        Check_DA{now,2}='日前市场发电机总出力上限拉格朗日乘子约束';
    end
end

for i=1:Num_Load                                        %日前市场用户每一时段每一用电段用电量互补松弛约束
    for j=1:Seg_bid
        for k=1:T
            now=now+1;
            Check_DA{now,1}=value(-gap<=(P_DA_D(i,j,k)-P_Dmin(i,j,k))<=M*c_DA_Dmin(i,j,k)+gap);
            Check_DA{now,2}='日前市场用户用电下限约束';
            
            now=now+1;
            Check_DA{now,1}=value(-gap<=u_DA_Dmin(i,j,k)<=M*(1-c_DA_Dmin(i,j,k))+gap);
            Check_DA{now,2}='日前市场用户用电下限拉格朗日乘子约束';
            
            now=now+1;
            Check_DA{now,1}=value(-gap<=(P_Dmax(i,j,k)-P_DA_D(i,j,k))<=M*c_DA_Dmax(i,j,k)+gap);
            Check_DA{now,2}='日前市场用户用电上限约束';
            
            now=now+1;
            Check_DA{now,1}=value(-gap<=u_DA_Dmax(i,j,k)<=M*(1-c_DA_Dmax(i,j,k))+gap);
            Check_DA{now,2}='日前市场用户用电上限拉格朗日乘子约束';
        end
    end
end

for i=1:Num_Node                                        %日前市场线路流过潮流互补松弛约束
    for j=1:size(AllNodes{i}.connectnodes,1)
        for k=1:T
            now=now+1;
            Check_DA{now,1}=value(-gap<=P_Lmax(i,AllNodes{i}.connectnodes(j))-B(i,AllNodes{i}.connectnodes(j))*(Delta_DA(i,k)-Delta_DA(AllNodes{i}.connectnodes(j),k))<=M*c_DA_Lmax(i,AllNodes{i}.connectnodes(j),k)+gap);
            Check_DA{now,2}='日前市场线路流过潮流约束';
            
            now=now+1;
            Check_DA{now,1}=value(-gap<=v_DA_Lmax(i,AllNodes{i}.connectnodes(j),k)<=M*(1-c_DA_Lmax(i,AllNodes{i}.connectnodes(j),k))+gap);
            Check_DA{now,2}='日前市场线路流过潮流拉格朗日乘子约束';
        end
    end
end

for i=1:Num_Node                                        %日前市场电压角互补松弛约束
    for j=1:T
        now=now+1;
        Check_DA{now,1}=value(-gap<=Delta_DA(i,j)+pi<=M*c_DA_Jmin(i,j)+gap);
        Check_DA{now,2}='日前市场电压角下限约束';
        
        now=now+1;
        Check_DA{now,1}=value(-gap<=e_DA_min(i,j)<=M*(1-c_DA_Jmin(i,j))+gap);
        Check_DA{now,2}='日前市场电压角下限拉格朗日乘子约束';
        
        now=now+1;
        Check_DA{now,1}=value(-gap<=pi-Delta_DA(i,j)<=M*c_DA_Jmax(i,j)+gap);
        Check_DA{now,2}='日前市场电压角上限约束';
        
        now=now+1;
        Check_DA{now,1}=value(-gap<=e_DA_max(i,j)<=M*(1-c_DA_Jmax(i,j))+gap);
        Check_DA{now,2}='日前市场电压角上限拉格朗日乘子约束';
    end
end

for i=1:Num_Node                                        %日前市场功率平衡约束
    for j=1:T
        New_Const=0;
        for k=1:Num_Gen
            if AllGens{k}.location==i
                New_Const=New_Const+sum(P_DA_G(k,:,j));
            end
        end
        for k=1:Num_Load
            if AllLoads{k}.location==i
                New_Const=New_Const-sum(P_DA_D(k,:,j));
            end
        end
        for k=1:size(AllNodes{i}.connectnodes,1)
            New_Const=New_Const-B(i,AllNodes{i}.connectnodes(k))*(Delta_DA(i,j)-Delta_DA(AllNodes{i}.connectnodes(k),j));
        end
        now=now+1;
        Check_DA{now,1}=(abs(value(New_Const))<=gap);
        Check_DA{now,2}='日前市场功率平衡约束';
    end
end

for i=1:T                                              %日前市场参考点电压角约束
    now=now+1;   
    Check_DA{now,1}=(abs(value(Delta_DA(1,i)))<=gap);
    Check_DA{now,2}='日前市场参考点电压角约束';
end

for i=1:Num_Const_DA
    if ~Check_DA{i,1}
        disp(['第',num2str(i),'条约束不满足要求，该约束的种类为',Check_DA{i,2}]);
        Iswrong_DA=1;
    end
end
if ~Iswrong_DA
    disp('日前市场约束符合要求');
end
%%
%实时市场约束正确性检验
Num_Const_RT=length(Const_RT);
Check_RT=cell(Num_Const_RT,2);
now=0; Iswrong_RT=0;
for m=1:W
    for i=1:Num_Gen                                          %KKT条件转换后约束1
        for j=1:Seg_offer
            for k=1:T
                now=now+1;
                if ismember(i,RE_Gens)
                    Check_RT{now,1}=(abs(value(Offer_Price(RE_Gens==i,j)-LMP_RT(AllGens{i}.location,m,k)+u_RT_Gmax(i,j,m,k)-u_RT_Gmin(i,j,m,k)+u_RT_GTmax(i,m,k)-u_RT_GTmin(i,m,k)))<=gap);
                    Check_RT{now,2}='KKT条件转换后约束1';
                else
                    Check_RT{now,1}=(abs(value(Marign_Cost(i,j)-LMP_RT(AllGens{i}.location,m,k)+u_RT_Gmax(i,j,m,k)-u_RT_Gmin(i,j,m,k)+u_RT_GTmax(i,m,k)-u_RT_GTmin(i,m,k)))<=gap);
                    Check_RT{now,2}='KKT条件转换后约束1';
                end
            end
        end
    end

    for i=1:Num_Load                                         %KKT条件转换后约束2
        for j=1:Seg_bid
            for k=1:T
                now=now+1;
                Check_RT{now,1}=(abs(value(LMP_RT(AllLoads{i}.location,m,k)-Bid_Price(i,j)+u_RT_Dmax(i,j,m,k)-u_RT_Dmin(i,j,m,k)))<=gap);
                Check_RT{now,2}='KKT条件转换后约束2';
            end
        end
    end

    for i=1:Num_Node                                         %KKT条件转换后约束3
        for j=1:T
            now=now+1;
            New_Const=0;
            for k=1:size(AllNodes{i}.connectnodes,1)
                New_Const=New_Const+B(i,AllNodes{i}.connectnodes(k))*(LMP_RT(i,m,j)-LMP_RT(AllNodes{i}.connectnodes(k),m,j))+B(i,AllNodes{i}.connectnodes(k))*(v_RT_Lmax(i,AllNodes{i}.connectnodes(k),m,j)-v_RT_Lmax(AllNodes{i}.connectnodes(k),i,m,j));
            end
            if i~=1
                Check_RT{now,1}=(abs(value(New_Const+e_RT_max(i,m,j)-e_RT_min(i,m,j)))<=gap);
                Check_RT{now,2}='KKT条件转换后约束3';
            else
                Check_RT{now,1}=(abs(value(New_Const+e_RT_max(i,m,j)-e_RT_min(i,m,j)+e_RT1(m,j)))<=gap);
                Check_RT{now,2}='KKT条件转换后约束3';
            end      
        end
    end

    for i=1:Num_Gen                                         %实时市场发电机每一时段每一报价段出力互补松弛约束
        for j=1:Seg_offer
            for k=1:T
                now=now+1;
                Check_RT{now,1}=value(-gap<=(P_RT_G(i,j,m,k)-P_RT_Gmin(i,j,m,k))<=M*c_RT_Gmin(i,j,m,k)+gap);
                Check_RT{now,2}='实时市场发电机出力下限约束';
                
                now=now+1;
                Check_RT{now,1}=value(-gap<=u_RT_Gmin(i,j,m,k)<=M*(1-c_RT_Gmin(i,j,m,k))+gap);
                Check_RT{now,2}='实时市场发电机出力下限拉格朗日乘子约束';
                
                now=now+1;
                Check_RT{now,1}=value(-gap<=(P_RT_Gmax(i,j,m,k)-P_RT_G(i,j,m,k))<=M*c_RT_Gmax(i,j,m,k)+gap);
                Check_RT{now,2}='实时市场发电机出力上限约束';
                
                now=now+1;
                Check_RT{now,1}=value(-gap<=u_RT_Gmax(i,j,m,k)<=M*(1-c_RT_Gmax(i,j,m,k))+gap);
                Check_RT{now,2}='实时市场发电机出力上限拉格朗日乘子约束';
            end
        end
    end

    for i=1:Num_Gen                                         %实时市场发电机每一时段出力互补松弛约束
        for j=1:T
            now=now+1;
            Check_RT{now,1}=value(-gap<=(sum(P_RT_G(i,:,m,j))-P_RT_GTmin(i,m,j))<=M*c_RT_GTmin(i,m,j)+gap);
            Check_RT{now,2}='实时市场发电机总出力下限约束';
            
            now=now+1;
            Check_RT{now,1}=value(-gap<=u_RT_GTmin(i,m,j)<=M*(1-c_RT_GTmin(i,m,j))+gap);
            Check_RT{now,2}='实时市场发电机总出力下限拉格朗日乘子约束';
            
            now=now+1;
            Check_RT{now,1}=value(-gap<=(P_RT_GTmax(i,m,j)-sum(P_RT_G(i,:,m,j)))<=M*c_RT_GTmax(i,m,j)+gap);
            Check_RT{now,2}='实时市场发电机总出力上限约束';
            
            now=now+1;
            Check_RT{now,1}=value(-gap<=u_RT_GTmax(i,m,j)<=M*(1-c_RT_GTmax(i,m,j))+gap);
            Check_RT{now,2}='实时市场发电机总出力上限拉格朗日乘子约束';
        end
    end

    for i=1:Num_Load                                        %实时市场用户每一时段每一用电段用电量互补松弛约束
        for j=1:Seg_bid
            for k=1:T
                now=now+1;
                Check_RT{now,1}=value(-gap<=(P_RT_D(i,j,m,k)-P_Dmin(i,j,k))<=M*c_RT_Dmin(i,j,m,k)+gap);
                Check_RT{now,2}='实时市场用户用电下限约束';
                
                now=now+1;
                Check_RT{now,1}=value(-gap<=u_RT_Dmin(i,j,m,k)<=M*(1-c_RT_Dmin(i,j,m,k))+gap);
                Check_RT{now,2}='实时市场用户用电下限拉格朗日乘子约束';
                
                now=now+1;
                Check_RT{now,1}=value(-gap<=(P_Dmax(i,j,k)-P_RT_D(i,j,m,k))<=M*c_RT_Dmax(i,j,m,k)+gap);
                Check_RT{now,2}='实时市场用户用电上限约束';
                
                now=now+1;
                Check_RT{now,1}=value(-gap<=u_RT_Dmax(i,j,m,k)<=M*(1-c_RT_Dmax(i,j,m,k))+gap);
                Check_RT{now,2}='实时市场用户用电上限拉格朗日乘子约束';
            end
        end
    end

    for i=1:Num_Node                                        %实时市场线路流过潮流互补松弛约束
        for j=1:size(AllNodes{i}.connectnodes,1)
            for k=1:T
                now=now+1;
                Check_RT{now,1}=value(-gap<=P_Lmax(i,AllNodes{i}.connectnodes(j))-B(i,AllNodes{i}.connectnodes(j))*(Delta_RT(i,m,k)-Delta_RT(AllNodes{i}.connectnodes(j),m,k))<=M*c_RT_Lmax(i,AllNodes{i}.connectnodes(j),m,k)+gap);
                Check_RT{now,2}='实时市场线路流过潮流约束';
                
                now=now+1;
                Check_RT{now,1}=value(-gap<=v_RT_Lmax(i,AllNodes{i}.connectnodes(j),m,k)<=M*(1-c_RT_Lmax(i,AllNodes{i}.connectnodes(j),m,k))+gap);
                Check_RT{now,2}='实时市场线路流过潮流拉格朗日乘子约束';
            end
        end
    end

    for i=1:Num_Node                                        %实时市场电压角互补松弛约束
        for j=1:T
            now=now+1;
            Check_RT{now,1}=value(-gap<=Delta_RT(i,m,j)+pi<=M*c_RT_Jmin(i,m,j)+gap);
            Check_RT{now,2}='实时市场电压角下限约束';
            
            now=now+1;
            Check_RT{now,1}=value(-gap<=e_RT_min(i,m,j)<=M*(1-c_RT_Jmin(i,m,j))+gap);
            Check_RT{now,2}='实时市场电压角下限拉格朗日乘子约束';
            
            now=now+1;
            Check_RT{now,1}=value(-gap<=pi-Delta_RT(i,m,j)<=M*c_RT_Jmax(i,m,j)+gap);
            Check_RT{now,2}='实时市场电压角上限约束';
            
            now=now+1;
            Check_RT{now,1}=value(-gap<=e_RT_max(i,m,j)<=M*(1-c_RT_Jmax(i,m,j))+gap);
            Check_RT{now,2}='实时市场电压角上限拉格朗日乘子约束';
        end
    end

    for i=1:Num_Node                                        %实时市场功率平衡约束
        for j=1:T
            New_Const=0;
            for k=1:Num_Gen
                if AllGens{k}.location==i
                    New_Const=New_Const+sum(P_RT_G(k,:,m,j));
                end
            end
            for k=1:Num_Load
                if AllLoads{k}.location==i
                    New_Const=New_Const-sum(P_RT_D(k,:,m,j));
                end
            end
            for k=1:size(AllNodes{i}.connectnodes,1)
                New_Const=New_Const-B(i,AllNodes{i}.connectnodes(k))*(Delta_RT(i,m,j)-Delta_RT(AllNodes{i}.connectnodes(k),m,j));
            end
            now=now+1;
            Check_RT{now,1}=(abs(value(New_Const))<=gap);
            Check_RT{now,2}='实时市场功率平衡约束';
        end
    end

    for i=1:T                                              %实时市场参考点电压角约束
        now=now+1;
        Check_RT{now,1}=(abs(value(Delta_RT(1,m,i)))<=gap);
        Check_RT{now,2}='实时市场参考点电压角约束';
    end
end

for i=1:Num_Const_RT
    if ~Check_RT{i,1}
        disp(['第',num2str(i),'条约束不满足要求，该约束的种类为',Check_RT{i,2}]);
        Iswrong_RT=1;
    end
end
if ~Iswrong_RT
    disp('实时市场约束符合要求');
end
