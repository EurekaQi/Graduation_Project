%%
clc,clear
cd 'D:\ѧϰ���ļ�\������\��Ŀ\��ҵ���\����&����'
%%
%��ȡ����
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
disp('���ݵ������');
%%
%������������
T=BasicPara{2,1};                             %����ʱ����
T_Gap=BasicPara{2,2};                         %����ÿ��ʱ����
W=BasicPara{2,17};                            %���ó�����
Seg_offer=BasicPara{2,3};                     %��������۶���
Seg_bid=BasicPara{2,4};                       %�û����۶���
Offer_UpperCup=BasicPara{2,7};                %�û���������
Offer_LowerCup=BasicPara{2,8};                %�û���������
Gama_FT=BasicPara{2,12};                      %δ��TGC�۸�����
Theta_FT=BasicPara{2,13};                     %TGC����ʱ�۸�����
Num_Approximate=BasicPara{2,14};              %���ƻ�׼ֵ����
RPS_ratio=BasicPara{2,15};                    %RPS�������
Penalty_Price=BasicPara{2,16};                %TGC����ĳͷ��۸�
RT_Wave=[0.992;1.025;1.012;0.994;1.001;1.017];%ʵʱ�г��еĲ���
RE_Gens=[7;8;10;12;14];                       %����Դ���繫˾��ӵ�еĻ�����
RE_Gens_Other=[9;11;13];                      %��������Դ��������
M=400000;
gap=1e-5;
Q_max=20;  Q_min=0;
Q_ref_gap=(Q_max-Q_min)/(Num_Approximate-1);                
%�ڵ���������
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
%�ڵ����ݴ���
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
%��·���ݴ���
B=zeros(size(NetBuses_Reform,1),size(NetBuses_Reform,1));
P_Lmax=zeros(size(NetBuses_Reform,1),size(NetBuses_Reform,1));
for i=1:size(NetBranches_Reform,1)
    B(NetBranches_Reform(i,1),NetBranches_Reform(i,2))=1/NetBranches_Reform(i,3);
    B(NetBranches_Reform(i,2),NetBranches_Reform(i,1))=1/NetBranches_Reform(i,3);
    P_Lmax(NetBranches_Reform(i,1),NetBranches_Reform(i,2))=NetBranches_Reform(i,4);
    P_Lmax(NetBranches_Reform(i,2),NetBranches_Reform(i,1))=NetBranches_Reform(i,4);
end
%�������ݴ���
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
%���������ݴ���
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
%��������������ݴ���
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
%�������������ݴ���
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
%�����л�������������һ��Ԫ��֮��
AllGens=[AllThermalGenerators;AllWindFarmGenerators;AllPhotovoltaicGenerators];   
%�������ݴ���
AllScenarios=cell(size(Scenarios,1)-1,1);
for i=1:size(Scenarios,1)-1
    probability=Scenarios{i+1,3};
    WF_change=Scenarios{i+1,4};
    PV_change=Scenarios{i+1,5};

    AllScenarios{i}=struct('probability',probability,'WF_change',WF_change,'PV_change',PV_change);
end
%�����������߼ʳɱ�ȷ��
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
%��ǰ�г������ÿ���۶γ���������ȷ��
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
%��ǰ�г������ÿһʱ��γ���������ȷ��
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
%ʵʱ�г������ÿ���۶γ���������ȷ��
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
%ʵʱ�г������ÿʱ��γ���������ȷ��
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
%�õ��û�ÿһ�α���������ȷ��
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
%�û��ı�������ȷ��
Bid_Price=zeros(size(AllLoads,1), Seg_bid);
for i=1:size(AllLoads,1)
    for j=1:Seg_bid
        Bid_Price(i,j)=BidData{i+1,j+1};
    end
end
%TGC�г�����ȷ��
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
%ʵʱ�г��������Ի�����ȷ��
P_max=zeros(size(RE_Gens,1),T); P_min=zeros(size(RE_Gens,1),T);
P_ref_gap=zeros(size(RE_Gens,1),T);
for i=1:size(RE_Gens,1)
    for j=1:T
        P_max(i,j)=P_DA_GTmax(RE_Gens(i),j);
        P_min(i,j)=P_DA_GTmin(RE_Gens(i),j);
        P_ref_gap(i,j)=(P_max(i,j)-P_min(i,j))/(Num_Approximate-1);
    end
end
disp('ԭʼ���ݴ������');
%%
%���һЩ��������
Num_Gen=size(AllGens,1);                                   %����ܷ������Ŀ���ڽ�ģ
Num_CoGen=size(RE_Gens,1);                                 %ս���Կ�������Դ���繫˾ӵ�еĻ�����Ŀ
Num_Load=size(AllLoads,1);                                 %����ܸ��������㽨ģ
Num_Node=size(AllNodes,1);                                 %�������·�����ڽ�ģ
%��������
Offer_Price=sdpvar(Num_CoGen,Seg_offer,'full');              %������������۱���
LMP_DA=sdpvar(Num_Node,T,'full');                          %��ǰ�г��еĽڵ���
P_DA_G=sdpvar(Num_Gen,Seg_offer,T,'full');                 %ÿһ̨���������ǰ�г���ͬʱ���Լ����ۿ�ĳ���
P_DA_D=sdpvar(Num_Load,Seg_bid,T,'full');                  %ÿһ���û�����ǰ�г���ͬʱ���Լ����ۿ���õ���
Delta_DA=sdpvar(Num_Node,T,'full');                        %��ǰ�г�ÿһ���ڵ�ĵ�ѹ��
LMP_RT=sdpvar(Num_Node,W,T,'full');                        %ʵʱ�г��еĽڵ���
P_RT_G=sdpvar(Num_Gen,Seg_offer,W,T,'full');               %ÿһ̨�������ʵʱ�г���ͬʱ�β�ͬ�����Լ����ۿ�ĳ���
P_RT_D=sdpvar(Num_Load,Seg_bid,W,T,'full');                %ÿһ���û���ʵʱ�г���ͬʱ�β�ͬ�����Լ����ۿ���õ���
Delta_RT=sdpvar(Num_Node,W,T,'full');                      %ʵʱ�г�ÿһ���ڵ�ĵ�ѹ��
TGC_Price_PT=sdpvar(W,1,'full');                           %��ǰ�г��е�TGC�۸�
TGC_Price_FT=sdpvar(1,1,'full');                           %δ���г��е�TGC�۸�
Q_PT=sdpvar(W,1,'full');                                   %�ڵ�ǰ�г��г��۵�TGC����
Q_PT_Other=sdpvar(W,1,'full');                             %������������Դ���繫˾�ڵ�ǰ�г����۵�TGC����
Q_FT=sdpvar(W,1,'full');                                   %��δ���г��г��۵�TGC����
Q_Gmax=sdpvar(W,1,'full');                                 %ս���Կ�������Դ���繫˾��ǰ�г��г��۵����TGC����
Q_Gmax_Other=sdpvar(W,1,'full');                           %������������Դ���繫˾��ǰ�г��г��۵����TGC����
%�������Ի������в����ı���
Z_DA1=sdpvar(T,1,'full');                                  %��ǰ�г��е����Ի���1
Z_DA2=sdpvar(T,1,'full');                                  %��ǰ�г��е����Ի���2
Z_RT11=sdpvar(W,T,'full');                                 %ʵʱ�г��е����Ի���1-1
Z_RT12=sdpvar(W,T,'full');                                 %ʵʱ�г��е����Ի���1-2
Z_RT2=sdpvar(W,T,'full');                                  %ʵʱ�г��е����Ի���2
Z_TGC=sdpvar(W,1,'full');                                  %TGC�г��е����Ի���
P_ref=sdpvar(Num_CoGen,Num_Approximate,'full');            %ʵʱ�г����Ի��еĽ��ƻ�׼ֵ
RT_Ancillary1=sdpvar(Num_CoGen,W,T,Num_Approximate,'full');%RT���Ի���������1
RT_Ancillary2=binvar(Num_CoGen,T,Num_Approximate,'full');  %RT���Ի���������2
Q_ref=sdpvar(Num_Approximate,1,'full');                    %TGC�г����Ի��еĽ��ƻ�׼ֵ
TGC_Ancillary1=sdpvar(W,Num_Approximate,'full');           %TGC���Ի���������1
TGC_Ancillary2=binvar(W,Num_Approximate,'full');           %TGC���Ի���������2
%����KKT�����е��������ճ��ӱ���
u_DA_Gmin=sdpvar(Num_Gen,Seg_offer,T,'full');              %��ǰ�г������ÿһʱ��ÿһ���۶γ���Լ������
u_DA_Gmax=sdpvar(Num_Gen,Seg_offer,T,'full');
u_DA_GTmin=sdpvar(Num_Gen,T,'full');                       %��ǰ�г������ÿһʱ�γ���Լ������
u_DA_GTmax=sdpvar(Num_Gen,T,'full');
u_DA_Dmin=sdpvar(Num_Load,Seg_bid,T,'full');               %��ǰ�г��û�ÿһʱ��ÿһ�õ���õ���Լ������
u_DA_Dmax=sdpvar(Num_Load,Seg_bid,T,'full');
v_DA_Lmax=sdpvar(Num_Node,Num_Node,T,'full');              %��ǰ�г���·��������Լ������
e_DA_min=sdpvar(Num_Node,T,'full');                        %��ǰ�г���ѹ��Լ������
e_DA_max=sdpvar(Num_Node,T,'full');
e_DA1=sdpvar(T,1,'full');                                  %��ǰ�г��ο��ڵ��ѹ��Լ������
u_RT_Gmin=sdpvar(Num_Gen,Seg_offer,W,T,'full');            %ʵʱ�г������ÿһʱ��ÿһ���۶�ÿһ��������Լ������
u_RT_Gmax=sdpvar(Num_Gen,Seg_offer,W,T,'full');
u_RT_GTmin=sdpvar(Num_Gen,W,T,'full');                     %ʵʱ�г������ÿһʱ��ÿһ��������Լ������
u_RT_GTmax=sdpvar(Num_Gen,W,T,'full');
u_RT_Dmin=sdpvar(Num_Load,Seg_bid,W,T,'full');             %ʵʱ�г��û�ÿһʱ��ÿһ�õ��ÿһ�����õ���Լ������
u_RT_Dmax=sdpvar(Num_Load,Seg_bid,W,T,'full');
v_RT_Lmax=sdpvar(Num_Node,Num_Node,W,T,'full');            %ʵʱ�г���·��������Լ������
e_RT_min=sdpvar(Num_Node,W,T,'full');                      %ʵʱ�г���ѹ��Լ������
e_RT_max=sdpvar(Num_Node,W,T,'full');
e_RT1=sdpvar(W,T,'full'); 
%����KKT�����ɳ�Լ�����Ի��еı���
c_DA_Gmin=binvar(Num_Gen,Seg_offer,T,'full');              %��ǰ�г������ÿһʱ��ÿһ���۶γ���Լ�����Ի�0-1����
c_DA_Gmax=binvar(Num_Gen,Seg_offer,T,'full');
c_DA_GTmin=binvar(Num_Gen,T,'full');                       %��ǰ�г������ÿһʱ�γ���Լ�����Ի�0-1����
c_DA_GTmax=binvar(Num_Gen,T,'full');
c_DA_Dmin=binvar(Num_Load,Seg_bid,T,'full');               %��ǰ�г��û�ÿһʱ��ÿһ�õ���õ���Լ�����Ի�0-1����
c_DA_Dmax=binvar(Num_Load,Seg_bid,T,'full');
c_DA_Lmax=binvar(Num_Node,Num_Node,T,'full');              %��ǰ�г���·��������Լ�����Ի�0-1����
c_DA_Jmin=binvar(Num_Node,T,'full');                       %��ǰ�г���ѹ��Լ�����Ի�0-1����
c_DA_Jmax=binvar(Num_Node,T,'full');
c_RT_Gmin=binvar(Num_Gen,Seg_offer,W,T,'full');            %ʵʱ�г������ÿһʱ��ÿһ���۶�ÿһ��������Լ�����Ի�0-1����
c_RT_Gmax=binvar(Num_Gen,Seg_offer,W,T,'full');
c_RT_GTmin=binvar(Num_Gen,W,T,'full');                     %ʵʱ�г������ÿһʱ��ÿһ��������Լ�����Ի�0-1����
c_RT_GTmax=binvar(Num_Gen,W,T,'full');
c_RT_Dmin=binvar(Num_Load,Seg_bid,W,T,'full');             %ʵʱ�г��û�ÿһʱ��ÿһ�õ��ÿһ�����õ���Լ�����Ի�0-1����
c_RT_Dmax=binvar(Num_Load,Seg_bid,W,T,'full');
c_RT_Lmax=binvar(Num_Node,Num_Node,W,T,'full');            %ʵʱ�г���·��������Լ�����Ի�0-1����
c_RT_Jmin=binvar(Num_Node,W,T,'full');                     %ʵʱ�г���ѹ��Լ�����Ի�0-1����
c_RT_Jmax=binvar(Num_Node,W,T,'full');
%����ɳ�
SC=sdpvar(1,1,'full');
disp('�����������');
%%
%����ʼֵ
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
%���Լ��
Const_Upper=[];
%��������Լ��     
for i=1:Num_CoGen                                         %����������Լ��            
    for j=1:Seg_offer
        Const_Upper=[Const_Upper,Offer_Price(i,j)<=Offer_UpperCup];
        Const_Upper=[Const_Upper,Offer_Price(i,j)>=Offer_LowerCup];
    end
end

for i=1:Num_CoGen                                         %���ۿ����Լ��            
    for j=2:Seg_offer
        Const_Upper=[Const_Upper,Offer_Price(i,j-1)<=Offer_Price(i,j)];
    end
end

for i=1:W                                                %TGC��������Լ��(������ʱ��Ϊ0)    
    Const_Upper=[Const_Upper,Q_PT(i)<=Q_Gmax(i)];
    Const_Upper=[Const_Upper,Q_PT(i)>=0];
end

for i=1:W                                                %ս���Կ�������Դ���繫˾TGC������������Լ��
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

for i=1:W                                                %������������Դ���繫˾TGC������������Լ��
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

for i=1:W                                                %Q_FT���㹫ʽԼ��
    Const_Upper=[Const_Upper,Q_FT(i)==Q_Gmax(i)-Q_PT(i)];
end

for i=1:W                                                %Q_PT_Other���㹫ʽԼ��
    Const_Upper=[Const_Upper,Q_PT_Other(i)==Q_Gmax_Other(i)];
end

disp('��������Լ��������');
%%
%�ײ�����һ����ǰ�г�����Լ��
Const_DA=[];

for i=1:Num_Gen                                          %KKT����ת����Լ��1
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

for i=1:Num_Load                                         %KKT����ת����Լ��2
    for j=1:Seg_bid
        for k=1:T
            Const_DA=[Const_DA,LMP_DA(AllLoads{i}.location,k)-Bid_Price(i,j)+u_DA_Dmax(i,j,k)-u_DA_Dmin(i,j,k)==0];
        end
    end
end

for i=1:Num_Node                                         %KKT����ת����Լ��3
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

for i=1:Num_Gen                                         %��ǰ�г������ÿһʱ��ÿһ���۶γ��������ɳ�Լ��
    for j=1:Seg_offer
        for k=1:T
            Const_DA=[Const_DA,0<=(P_DA_G(i,j,k)-P_DA_Gmin(i,j,k))<=M*c_DA_Gmin(i,j,k)];
            Const_DA=[Const_DA,0<=u_DA_Gmin(i,j,k)<=M*(1-c_DA_Gmin(i,j,k))];
            Const_DA=[Const_DA,0<=(P_DA_Gmax(i,j,k)-P_DA_G(i,j,k))<=M*c_DA_Gmax(i,j,k)];
            Const_DA=[Const_DA,0<=u_DA_Gmax(i,j,k)<=M*(1-c_DA_Gmax(i,j,k))];
        end
    end
end

for i=1:Num_Gen                                         %��ǰ�г������ÿһʱ�γ��������ɳ�Լ��
    for j=1:T
        Const_DA=[Const_DA,0<=(sum(P_DA_G(i,:,j))-P_DA_GTmin(i,j))<=M*c_DA_GTmin(i,j)];
        Const_DA=[Const_DA,0<=u_DA_GTmin(i,j)<=M*(1-c_DA_GTmin(i,j))];
        Const_DA=[Const_DA,0<=(P_DA_GTmax(i,j)-sum(P_DA_G(i,:,j)))<=M*c_DA_GTmax(i,j)];
        Const_DA=[Const_DA,0<=u_DA_GTmax(i,j)<=M*(1-c_DA_GTmax(i,j))];
    end
end

for i=1:Num_Load                                        %��ǰ�г��û�ÿһʱ��ÿһ�õ���õ��������ɳ�Լ��
    for j=1:Seg_bid
        for k=1:T
            Const_DA=[Const_DA,0<=(P_DA_D(i,j,k)-P_Dmin(i,j,k))<=M*c_DA_Dmin(i,j,k)];
            Const_DA=[Const_DA,0<=u_DA_Dmin(i,j,k)<=M*(1-c_DA_Dmin(i,j,k))];
            Const_DA=[Const_DA,0<=(P_Dmax(i,j,k)-P_DA_D(i,j,k))<=M*c_DA_Dmax(i,j,k)];
            Const_DA=[Const_DA,0<=u_DA_Dmax(i,j,k)<=M*(1-c_DA_Dmax(i,j,k))];
        end
    end
end

for i=1:Num_Node                                        %��ǰ�г���·�������������ɳ�Լ��
    for j=1:size(AllNodes{i}.connectnodes,1)
        for k=1:T
            Const_DA=[Const_DA,0<=P_Lmax(i,AllNodes{i}.connectnodes(j))-B(i,AllNodes{i}.connectnodes(j))*(Delta_DA(i,k)-Delta_DA(AllNodes{i}.connectnodes(j),k))<=M*c_DA_Lmax(i,AllNodes{i}.connectnodes(j),k)];
            Const_DA=[Const_DA,0<=v_DA_Lmax(i,AllNodes{i}.connectnodes(j),k)<=M*(1-c_DA_Lmax(i,AllNodes{i}.connectnodes(j),k))];
        end
    end
end

for i=1:Num_Node                                        %��ǰ�г���ѹ�ǻ����ɳ�Լ��
    for j=1:T
        Const_DA=[Const_DA,0<=Delta_DA(i,j)+pi<=M*c_DA_Jmin(i,j)];
        Const_DA=[Const_DA,0<=e_DA_min(i,j)<=M*(1-c_DA_Jmin(i,j))];
        Const_DA=[Const_DA,0<=pi-Delta_DA(i,j)<=M*c_DA_Jmax(i,j)];
        Const_DA=[Const_DA,0<=e_DA_max(i,j)<=M*(1-c_DA_Jmax(i,j))];
    end
end

for i=1:Num_Node                                        %��ǰ�г�����ƽ��Լ��
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

for i=1:T                                              %��ǰ�г��ο����ѹ��Լ��
    Const_DA=[Const_DA,Delta_DA(1,i)==0];
end
disp('�ײ�����һ����Լ��������');
%%
%�ײ��������ʵʱ�г�����Լ��
Const_RT=[];
for m=1:W
    for i=1:Num_Gen                                          %KKT����ת����Լ��1
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

    for i=1:Num_Load                                         %KKT����ת����Լ��2
        for j=1:Seg_bid
            for k=1:T
                Const_RT=[Const_RT,LMP_RT(AllLoads{i}.location,m,k)-Bid_Price(i,j)+u_RT_Dmax(i,j,m,k)-u_RT_Dmin(i,j,m,k)==0];
            end
        end
    end

    for i=1:Num_Node                                         %KKT����ת����Լ��3
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

    for i=1:Num_Gen                                         %ʵʱ�г������ÿһʱ��ÿһ���۶γ��������ɳ�Լ��
        for j=1:Seg_offer
            for k=1:T
                Const_RT=[Const_RT,0<=(P_RT_G(i,j,m,k)-P_RT_Gmin(i,j,m,k))<=M*c_RT_Gmin(i,j,m,k)];
                Const_RT=[Const_RT,0<=u_RT_Gmin(i,j,m,k)<=M*(1-c_RT_Gmin(i,j,m,k))];
                Const_RT=[Const_RT,0<=(P_RT_Gmax(i,j,m,k)-P_RT_G(i,j,m,k))<=M*c_RT_Gmax(i,j,m,k)];
                Const_RT=[Const_RT,0<=u_RT_Gmax(i,j,m,k)<=M*(1-c_RT_Gmax(i,j,m,k))];
            end
        end
    end

    for i=1:Num_Gen                                         %ʵʱ�г������ÿһʱ�γ��������ɳ�Լ��
        for j=1:T
            Const_RT=[Const_RT,0<=(sum(P_RT_G(i,:,m,j))-P_RT_GTmin(i,m,j))<=M*c_RT_GTmin(i,m,j)];
            Const_RT=[Const_RT,0<=u_RT_GTmin(i,m,j)<=M*(1-c_RT_GTmin(i,m,j))];
            Const_RT=[Const_RT,0<=(P_RT_GTmax(i,m,j)-sum(P_RT_G(i,:,m,j)))<=M*c_RT_GTmax(i,m,j)];
            Const_RT=[Const_RT,0<=u_RT_GTmax(i,m,j)<=M*(1-c_RT_GTmax(i,m,j))];
        end
    end

    for i=1:Num_Load                                        %ʵʱ�г��û�ÿһʱ��ÿһ�õ���õ��������ɳ�Լ��
        for j=1:Seg_bid
            for k=1:T
                Const_RT=[Const_RT,0<=(P_RT_D(i,j,m,k)-P_Dmin(i,j,k))<=M*c_RT_Dmin(i,j,m,k)];
                Const_RT=[Const_RT,0<=u_RT_Dmin(i,j,m,k)<=M*(1-c_RT_Dmin(i,j,m,k))];
                Const_RT=[Const_RT,0<=(P_Dmax(i,j,k)-P_RT_D(i,j,m,k))<=M*c_RT_Dmax(i,j,m,k)];
                Const_RT=[Const_RT,0<=u_RT_Dmax(i,j,m,k)<=M*(1-c_RT_Dmax(i,j,m,k))];
            end
        end
    end

    for i=1:Num_Node                                        %ʵʱ�г���·�������������ɳ�Լ��
        for j=1:size(AllNodes{i}.connectnodes,1)
            for k=1:T
                Const_RT=[Const_RT,0<=P_Lmax(i,AllNodes{i}.connectnodes(j))-B(i,AllNodes{i}.connectnodes(j))*(Delta_RT(i,m,k)-Delta_RT(AllNodes{i}.connectnodes(j),m,k))<=M*c_RT_Lmax(i,AllNodes{i}.connectnodes(j),m,k)];
                Const_RT=[Const_RT,0<=v_RT_Lmax(i,AllNodes{i}.connectnodes(j),m,k)<=M*(1-c_RT_Lmax(i,AllNodes{i}.connectnodes(j),m,k))];
            end
        end
    end

    for i=1:Num_Node                                        %ʵʱ�г���ѹ�ǻ����ɳ�Լ��
        for j=1:T
            Const_RT=[Const_RT,0<=Delta_RT(i,m,j)+pi<=M*c_RT_Jmin(i,m,j)];
            Const_RT=[Const_RT,0<=e_RT_min(i,m,j)<=M*(1-c_RT_Jmin(i,m,j))];
            Const_RT=[Const_RT,0<=pi-Delta_RT(i,m,j)<=M*c_RT_Jmax(i,m,j)];
            Const_RT=[Const_RT,0<=e_RT_max(i,m,j)<=M*(1-c_RT_Jmax(i,m,j))];
        end
    end

    for i=1:Num_Node                                        %ʵʱ�г�����ƽ��Լ��
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

    for i=1:T                                              %ʵʱ�г��ο����ѹ��Լ��
        Const_RT=[Const_RT,Delta_RT(1,m,i)==0];
    end
end
disp('�ײ����������Լ��������');
%%
Const_Liner_DA=[];

%��ǰ�г��������Ի�Լ��
for i=1:T                                            %��ǰ�г����Ի���Z_DA1���㹫ʽ
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

for i=1:T                                            %��ǰ�г����Ի���Z_DA2���㹫ʽ
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
disp('��ǰ�г����Ի���Լ��������');
%%
Const_Liner_RT=[];

%ʵʱ�г��������Ի�Լ��
for i=1:W                                            %ʵʱ�г����Ի���Z_RT11���㹫ʽ
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

for i=1:W                                            %ʵʱ�г����Ի���Z_RT12���㹫ʽ        
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

for i=1:W                                            %ʵʱ�г����Ի���Z_RT2���㹫ʽ        
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
disp('ʵʱ�г����Ի���Լ��������');
%%
Const_Liner_TGC=[];

%TGC�г��������Ի�Լ��
for i=1:W                                           %TGC�г�Z_TGC���㹫ʽ
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
disp('TGC�г����Ի���Լ��������');
%%
%����Ŀ�꺯��
R_DA=(sum(Z_DA1)+sum(Z_DA2))*T_Gap;                       %��ǰ�г�������� 

New_Const=0;                                              %ʵʱ�г��������      
for i=1:W
    New_Const=New_Const+AllScenarios{i}.probability*(sum(Z_RT11(i,:))+sum(Z_RT12(i,:))+sum(Z_RT2(i,:)))*T_Gap;   
end
R_RT=New_Const;

New_Const=0;                                              %TGC�г��������Լ��       
for i=1:W
    New_Const=New_Const+AllScenarios{i}.probability*(alpha_TGC*Q_PT(i)-beta_TGC*Z_TGC(i)-beta_TGC*Q_PT(i)*Q_PT(i)+Gama_FT*Penalty_Price*(Q_Gmax(i)-Q_PT(i)));   
end
R_TGC=New_Const;

Obj=-(R_DA+R_RT+R_TGC);
%��ʼ���
Ops=sdpsettings('verbose',1,'solver','gurobi','usex0',1);    %����������
%Const=[Const_Upper,Const_DA,Const_RT,Const_Liner_DA,Const_Liner_RT,Const_Liner_TGC];
disp('��ʼ���');

Result=solvesdp([Const_Upper,Const_RT],0,Ops);
if Result.problem==0 
    fprintf('���ɹ��������ϢΪ��%s\n',Result.info);
else
    fprintf('�������г���������ϢΪ��%s\n',Result.info);
end
%%
%��ǰ�г�Լ����ȷ�Լ���
Num_Const_DA=length(Const_DA);
Check_DA=cell(Num_Const_DA,2);
now=0; Iswrong_DA=0;
for i=1:Num_Gen                                          %KKT����ת����Լ��1
    for j=1:Seg_offer
        for k=1:T
            now=now+1;
            if ismember(i,RE_Gens)
                Check_DA{now,1}=(abs(value(Offer_Price(RE_Gens==i,j)-LMP_DA(AllGens{i}.location,k)+u_DA_Gmax(i,j,k)-u_DA_Gmin(i,j,k)+u_DA_GTmax(i,k)-u_DA_GTmin(i,k)))<=gap);
                Check_DA{now,2}='KKT����ת����Լ��1';
            else
                Check_DA{now,1}=(abs(value(Marign_Cost(i,j)-LMP_DA(AllGens{i}.location,k)+u_DA_Gmax(i,j,k)-u_DA_Gmin(i,j,k)+u_DA_GTmax(i,k)-u_DA_GTmin(i,k)))<=gap);
                Check_DA{now,2}='KKT����ת����Լ��1';
            end
        end
    end
end

for i=1:Num_Load                                         %KKT����ת����Լ��2
    for j=1:Seg_bid
        for k=1:T
            now=now+1;
            Check_DA{now,1}=(abs(value(LMP_DA(AllLoads{i}.location,k)-Bid_Price(i,j)+u_DA_Dmax(i,j,k)-u_DA_Dmin(i,j,k)))<=gap);
            Check_DA{now,2}='KKT����ת����Լ��2';
        end
    end
end

for i=1:Num_Node                                         %KKT����ת����Լ��3
    for j=1:T
        now=now+1;
        New_Const=0;
        for k=1:size(AllNodes{i}.connectnodes,1)
            New_Const=New_Const+B(i,AllNodes{i}.connectnodes(k))*(LMP_DA(i,j)-LMP_DA(AllNodes{i}.connectnodes(k),j))+B(i,AllNodes{i}.connectnodes(k))*(v_DA_Lmax(i,AllNodes{i}.connectnodes(k),j)-v_DA_Lmax(AllNodes{i}.connectnodes(k),i,j));
        end
        if i~=1
            Check_DA{now,1}=(abs(value(New_Const+e_DA_max(i,j)-e_DA_min(i,j)))<=gap);
            Check_DA{now,2}='KKT����ת����Լ��3';
        else
            Check_DA{now,1}=(abs(value(New_Const+e_DA_max(i,j)-e_DA_min(i,j)+e_DA1(j)))<=gap);
            Check_DA{now,2}='KKT����ת����Լ��3';
        end      
    end
end

for i=1:Num_Gen                                         %��ǰ�г������ÿһʱ��ÿһ���۶γ��������ɳ�Լ��
    for j=1:Seg_offer
        for k=1:T
            now=now+1;
            Check_DA{now,1}=value(-gap<=(P_DA_G(i,j,k)-P_DA_Gmin(i,j,k))<=M*c_DA_Gmin(i,j,k)+gap);
            Check_DA{now,2}='��ǰ�г��������������Լ��';
            
            now=now+1;
            Check_DA{now,1}=value(-gap<=u_DA_Gmin(i,j,k)<=M*(1-c_DA_Gmin(i,j,k))+gap);
            Check_DA{now,2}='��ǰ�г���������������������ճ���Լ��';
            
            now=now+1;
            Check_DA{now,1}=value(-gap<=(P_DA_Gmax(i,j,k)-P_DA_G(i,j,k))<=M*c_DA_Gmax(i,j,k)+gap);
            Check_DA{now,2}='��ǰ�г��������������Լ��';
            
            now=now+1;
            Check_DA{now,1}=value(-gap<=u_DA_Gmax(i,j,k)<=M*(1-c_DA_Gmax(i,j,k))+gap);
            Check_DA{now,2}='��ǰ�г���������������������ճ���Լ��';
        end
    end
end

for i=1:Num_Gen                                         %��ǰ�г������ÿһʱ�γ��������ɳ�Լ��
    for j=1:T
        now=now+1;
        Check_DA{now,1}=value(-gap<=(sum(P_DA_G(i,:,j))-P_DA_GTmin(i,j))<=M*c_DA_GTmin(i,j)+gap);
        Check_DA{now,2}='��ǰ�г�������ܳ�������Լ��';
        
        now=now+1;
        Check_DA{now,1}=value(-gap<=u_DA_GTmin(i,j)<=M*(1-c_DA_GTmin(i,j))+gap);
        Check_DA{now,2}='��ǰ�г�������ܳ��������������ճ���Լ��';
        
        now=now+1;
        Check_DA{now,1}=value(-gap<=(P_DA_GTmax(i,j)-sum(P_DA_G(i,:,j)))<=M*c_DA_GTmax(i,j)+gap);
        Check_DA{now,2}='��ǰ�г�������ܳ�������Լ��';
        
        now=now+1;
        Check_DA{now,1}=value(-gap<=u_DA_GTmax(i,j)<=M*(1-c_DA_GTmax(i,j))+gap);
        Check_DA{now,2}='��ǰ�г�������ܳ��������������ճ���Լ��';
    end
end

for i=1:Num_Load                                        %��ǰ�г��û�ÿһʱ��ÿһ�õ���õ��������ɳ�Լ��
    for j=1:Seg_bid
        for k=1:T
            now=now+1;
            Check_DA{now,1}=value(-gap<=(P_DA_D(i,j,k)-P_Dmin(i,j,k))<=M*c_DA_Dmin(i,j,k)+gap);
            Check_DA{now,2}='��ǰ�г��û��õ�����Լ��';
            
            now=now+1;
            Check_DA{now,1}=value(-gap<=u_DA_Dmin(i,j,k)<=M*(1-c_DA_Dmin(i,j,k))+gap);
            Check_DA{now,2}='��ǰ�г��û��õ������������ճ���Լ��';
            
            now=now+1;
            Check_DA{now,1}=value(-gap<=(P_Dmax(i,j,k)-P_DA_D(i,j,k))<=M*c_DA_Dmax(i,j,k)+gap);
            Check_DA{now,2}='��ǰ�г��û��õ�����Լ��';
            
            now=now+1;
            Check_DA{now,1}=value(-gap<=u_DA_Dmax(i,j,k)<=M*(1-c_DA_Dmax(i,j,k))+gap);
            Check_DA{now,2}='��ǰ�г��û��õ������������ճ���Լ��';
        end
    end
end

for i=1:Num_Node                                        %��ǰ�г���·�������������ɳ�Լ��
    for j=1:size(AllNodes{i}.connectnodes,1)
        for k=1:T
            now=now+1;
            Check_DA{now,1}=value(-gap<=P_Lmax(i,AllNodes{i}.connectnodes(j))-B(i,AllNodes{i}.connectnodes(j))*(Delta_DA(i,k)-Delta_DA(AllNodes{i}.connectnodes(j),k))<=M*c_DA_Lmax(i,AllNodes{i}.connectnodes(j),k)+gap);
            Check_DA{now,2}='��ǰ�г���·��������Լ��';
            
            now=now+1;
            Check_DA{now,1}=value(-gap<=v_DA_Lmax(i,AllNodes{i}.connectnodes(j),k)<=M*(1-c_DA_Lmax(i,AllNodes{i}.connectnodes(j),k))+gap);
            Check_DA{now,2}='��ǰ�г���·���������������ճ���Լ��';
        end
    end
end

for i=1:Num_Node                                        %��ǰ�г���ѹ�ǻ����ɳ�Լ��
    for j=1:T
        now=now+1;
        Check_DA{now,1}=value(-gap<=Delta_DA(i,j)+pi<=M*c_DA_Jmin(i,j)+gap);
        Check_DA{now,2}='��ǰ�г���ѹ������Լ��';
        
        now=now+1;
        Check_DA{now,1}=value(-gap<=e_DA_min(i,j)<=M*(1-c_DA_Jmin(i,j))+gap);
        Check_DA{now,2}='��ǰ�г���ѹ�������������ճ���Լ��';
        
        now=now+1;
        Check_DA{now,1}=value(-gap<=pi-Delta_DA(i,j)<=M*c_DA_Jmax(i,j)+gap);
        Check_DA{now,2}='��ǰ�г���ѹ������Լ��';
        
        now=now+1;
        Check_DA{now,1}=value(-gap<=e_DA_max(i,j)<=M*(1-c_DA_Jmax(i,j))+gap);
        Check_DA{now,2}='��ǰ�г���ѹ�������������ճ���Լ��';
    end
end

for i=1:Num_Node                                        %��ǰ�г�����ƽ��Լ��
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
        Check_DA{now,2}='��ǰ�г�����ƽ��Լ��';
    end
end

for i=1:T                                              %��ǰ�г��ο����ѹ��Լ��
    now=now+1;   
    Check_DA{now,1}=(abs(value(Delta_DA(1,i)))<=gap);
    Check_DA{now,2}='��ǰ�г��ο����ѹ��Լ��';
end

for i=1:Num_Const_DA
    if ~Check_DA{i,1}
        disp(['��',num2str(i),'��Լ��������Ҫ�󣬸�Լ��������Ϊ',Check_DA{i,2}]);
        Iswrong_DA=1;
    end
end
if ~Iswrong_DA
    disp('��ǰ�г�Լ������Ҫ��');
end
%%
%ʵʱ�г�Լ����ȷ�Լ���
Num_Const_RT=length(Const_RT);
Check_RT=cell(Num_Const_RT,2);
now=0; Iswrong_RT=0;
for m=1:W
    for i=1:Num_Gen                                          %KKT����ת����Լ��1
        for j=1:Seg_offer
            for k=1:T
                now=now+1;
                if ismember(i,RE_Gens)
                    Check_RT{now,1}=(abs(value(Offer_Price(RE_Gens==i,j)-LMP_RT(AllGens{i}.location,m,k)+u_RT_Gmax(i,j,m,k)-u_RT_Gmin(i,j,m,k)+u_RT_GTmax(i,m,k)-u_RT_GTmin(i,m,k)))<=gap);
                    Check_RT{now,2}='KKT����ת����Լ��1';
                else
                    Check_RT{now,1}=(abs(value(Marign_Cost(i,j)-LMP_RT(AllGens{i}.location,m,k)+u_RT_Gmax(i,j,m,k)-u_RT_Gmin(i,j,m,k)+u_RT_GTmax(i,m,k)-u_RT_GTmin(i,m,k)))<=gap);
                    Check_RT{now,2}='KKT����ת����Լ��1';
                end
            end
        end
    end

    for i=1:Num_Load                                         %KKT����ת����Լ��2
        for j=1:Seg_bid
            for k=1:T
                now=now+1;
                Check_RT{now,1}=(abs(value(LMP_RT(AllLoads{i}.location,m,k)-Bid_Price(i,j)+u_RT_Dmax(i,j,m,k)-u_RT_Dmin(i,j,m,k)))<=gap);
                Check_RT{now,2}='KKT����ת����Լ��2';
            end
        end
    end

    for i=1:Num_Node                                         %KKT����ת����Լ��3
        for j=1:T
            now=now+1;
            New_Const=0;
            for k=1:size(AllNodes{i}.connectnodes,1)
                New_Const=New_Const+B(i,AllNodes{i}.connectnodes(k))*(LMP_RT(i,m,j)-LMP_RT(AllNodes{i}.connectnodes(k),m,j))+B(i,AllNodes{i}.connectnodes(k))*(v_RT_Lmax(i,AllNodes{i}.connectnodes(k),m,j)-v_RT_Lmax(AllNodes{i}.connectnodes(k),i,m,j));
            end
            if i~=1
                Check_RT{now,1}=(abs(value(New_Const+e_RT_max(i,m,j)-e_RT_min(i,m,j)))<=gap);
                Check_RT{now,2}='KKT����ת����Լ��3';
            else
                Check_RT{now,1}=(abs(value(New_Const+e_RT_max(i,m,j)-e_RT_min(i,m,j)+e_RT1(m,j)))<=gap);
                Check_RT{now,2}='KKT����ת����Լ��3';
            end      
        end
    end

    for i=1:Num_Gen                                         %ʵʱ�г������ÿһʱ��ÿһ���۶γ��������ɳ�Լ��
        for j=1:Seg_offer
            for k=1:T
                now=now+1;
                Check_RT{now,1}=value(-gap<=(P_RT_G(i,j,m,k)-P_RT_Gmin(i,j,m,k))<=M*c_RT_Gmin(i,j,m,k)+gap);
                Check_RT{now,2}='ʵʱ�г��������������Լ��';
                
                now=now+1;
                Check_RT{now,1}=value(-gap<=u_RT_Gmin(i,j,m,k)<=M*(1-c_RT_Gmin(i,j,m,k))+gap);
                Check_RT{now,2}='ʵʱ�г���������������������ճ���Լ��';
                
                now=now+1;
                Check_RT{now,1}=value(-gap<=(P_RT_Gmax(i,j,m,k)-P_RT_G(i,j,m,k))<=M*c_RT_Gmax(i,j,m,k)+gap);
                Check_RT{now,2}='ʵʱ�г��������������Լ��';
                
                now=now+1;
                Check_RT{now,1}=value(-gap<=u_RT_Gmax(i,j,m,k)<=M*(1-c_RT_Gmax(i,j,m,k))+gap);
                Check_RT{now,2}='ʵʱ�г���������������������ճ���Լ��';
            end
        end
    end

    for i=1:Num_Gen                                         %ʵʱ�г������ÿһʱ�γ��������ɳ�Լ��
        for j=1:T
            now=now+1;
            Check_RT{now,1}=value(-gap<=(sum(P_RT_G(i,:,m,j))-P_RT_GTmin(i,m,j))<=M*c_RT_GTmin(i,m,j)+gap);
            Check_RT{now,2}='ʵʱ�г�������ܳ�������Լ��';
            
            now=now+1;
            Check_RT{now,1}=value(-gap<=u_RT_GTmin(i,m,j)<=M*(1-c_RT_GTmin(i,m,j))+gap);
            Check_RT{now,2}='ʵʱ�г�������ܳ��������������ճ���Լ��';
            
            now=now+1;
            Check_RT{now,1}=value(-gap<=(P_RT_GTmax(i,m,j)-sum(P_RT_G(i,:,m,j)))<=M*c_RT_GTmax(i,m,j)+gap);
            Check_RT{now,2}='ʵʱ�г�������ܳ�������Լ��';
            
            now=now+1;
            Check_RT{now,1}=value(-gap<=u_RT_GTmax(i,m,j)<=M*(1-c_RT_GTmax(i,m,j))+gap);
            Check_RT{now,2}='ʵʱ�г�������ܳ��������������ճ���Լ��';
        end
    end

    for i=1:Num_Load                                        %ʵʱ�г��û�ÿһʱ��ÿһ�õ���õ��������ɳ�Լ��
        for j=1:Seg_bid
            for k=1:T
                now=now+1;
                Check_RT{now,1}=value(-gap<=(P_RT_D(i,j,m,k)-P_Dmin(i,j,k))<=M*c_RT_Dmin(i,j,m,k)+gap);
                Check_RT{now,2}='ʵʱ�г��û��õ�����Լ��';
                
                now=now+1;
                Check_RT{now,1}=value(-gap<=u_RT_Dmin(i,j,m,k)<=M*(1-c_RT_Dmin(i,j,m,k))+gap);
                Check_RT{now,2}='ʵʱ�г��û��õ������������ճ���Լ��';
                
                now=now+1;
                Check_RT{now,1}=value(-gap<=(P_Dmax(i,j,k)-P_RT_D(i,j,m,k))<=M*c_RT_Dmax(i,j,m,k)+gap);
                Check_RT{now,2}='ʵʱ�г��û��õ�����Լ��';
                
                now=now+1;
                Check_RT{now,1}=value(-gap<=u_RT_Dmax(i,j,m,k)<=M*(1-c_RT_Dmax(i,j,m,k))+gap);
                Check_RT{now,2}='ʵʱ�г��û��õ������������ճ���Լ��';
            end
        end
    end

    for i=1:Num_Node                                        %ʵʱ�г���·�������������ɳ�Լ��
        for j=1:size(AllNodes{i}.connectnodes,1)
            for k=1:T
                now=now+1;
                Check_RT{now,1}=value(-gap<=P_Lmax(i,AllNodes{i}.connectnodes(j))-B(i,AllNodes{i}.connectnodes(j))*(Delta_RT(i,m,k)-Delta_RT(AllNodes{i}.connectnodes(j),m,k))<=M*c_RT_Lmax(i,AllNodes{i}.connectnodes(j),m,k)+gap);
                Check_RT{now,2}='ʵʱ�г���·��������Լ��';
                
                now=now+1;
                Check_RT{now,1}=value(-gap<=v_RT_Lmax(i,AllNodes{i}.connectnodes(j),m,k)<=M*(1-c_RT_Lmax(i,AllNodes{i}.connectnodes(j),m,k))+gap);
                Check_RT{now,2}='ʵʱ�г���·���������������ճ���Լ��';
            end
        end
    end

    for i=1:Num_Node                                        %ʵʱ�г���ѹ�ǻ����ɳ�Լ��
        for j=1:T
            now=now+1;
            Check_RT{now,1}=value(-gap<=Delta_RT(i,m,j)+pi<=M*c_RT_Jmin(i,m,j)+gap);
            Check_RT{now,2}='ʵʱ�г���ѹ������Լ��';
            
            now=now+1;
            Check_RT{now,1}=value(-gap<=e_RT_min(i,m,j)<=M*(1-c_RT_Jmin(i,m,j))+gap);
            Check_RT{now,2}='ʵʱ�г���ѹ�������������ճ���Լ��';
            
            now=now+1;
            Check_RT{now,1}=value(-gap<=pi-Delta_RT(i,m,j)<=M*c_RT_Jmax(i,m,j)+gap);
            Check_RT{now,2}='ʵʱ�г���ѹ������Լ��';
            
            now=now+1;
            Check_RT{now,1}=value(-gap<=e_RT_max(i,m,j)<=M*(1-c_RT_Jmax(i,m,j))+gap);
            Check_RT{now,2}='ʵʱ�г���ѹ�������������ճ���Լ��';
        end
    end

    for i=1:Num_Node                                        %ʵʱ�г�����ƽ��Լ��
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
            Check_RT{now,2}='ʵʱ�г�����ƽ��Լ��';
        end
    end

    for i=1:T                                              %ʵʱ�г��ο����ѹ��Լ��
        now=now+1;
        Check_RT{now,1}=(abs(value(Delta_RT(1,m,i)))<=gap);
        Check_RT{now,2}='ʵʱ�г��ο����ѹ��Լ��';
    end
end

for i=1:Num_Const_RT
    if ~Check_RT{i,1}
        disp(['��',num2str(i),'��Լ��������Ҫ�󣬸�Լ��������Ϊ',Check_RT{i,2}]);
        Iswrong_RT=1;
    end
end
if ~Iswrong_RT
    disp('ʵʱ�г�Լ������Ҫ��');
end
