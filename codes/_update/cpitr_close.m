
% CPITR_close.m

%%Analyze stability in the close case for the (current 
%and expected) CPI inflation Taylor Rule

clear;
clear all 
global sigma shi chi fi beta eta gama eps sc tita v delta landa k1 k2 w mu si

%********* Initialize fundamental Parameters **************************************
sigma=5;        %IES consumption {1,5}
shi=-1;         %Share desutility labor
chi=0.47;       %Frisch elasticity {0.47,3}
fi=10^(-6);  %Endogenerity of discounting {0,10^(-6)}
beta=0.99;      %stationary discounting
raro=0;         %parameter to calibrate discounting
eta=1.5;        %Elasticity substitution in composite consumption {1,1.5}
gama=0;       %share in composite consumption (openess in demand) {0,0.4}
eps=6;          %Elasticity substitution in agregator consumption
sc=0.8;         %Consumption to output ratio
tita=0.75;      %Calvo probability
v=0;           %Elasticity of substititution in production {0,-2}
delta=0;    %Openess in production {0,0.144}

%*********Reduced parameter**********
tau=sigma*(1-gama)/((1-gama)^2+sigma*gama*eta*(2-gama));
%New keynesian Phillips
landa=((1-tita)*(1-tita*beta)/tita)*((1-v)*(1-delta)/(1-v+...
    delta*chi));
k1=(sigma/(1-gama))+chi;
%control to reduce to complete market
k2=0;
%Is equation
w=sigma/(sigma-fi);
mu=((1-gama)/(sigma-fi))*(1-gama+eta*gama*(2-gama)*(sigma-fi)/(1-gama));
si=eta*gama*fi*(2-gama)/((1-gama)*(sigma-fi));

%******THREE CASES TO COMPARE******
%we will analyze the imcomplete case, which is the baseline, and then we
%compare with other two cases: (I) Comlete case (kdos=0), (II) and the
%close economy (gama=0, and delta=0).

%**********Stability analysis***********
decimals=101;    
fixmax=4;
fipimax=4;
fis=0;
m=fixmax*decimals;  %interaciones en fix 
n=fipimax*decimals; %interaciones en fipi 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%:::::::::::CPI Inflation Taylor Rule (CPITR):::::::::::::::::::::::
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creating the output matrix
fixs=zeros(m,1); %parameter of output gap
fipis=zeros(n,1);   %parameter of inflation gap
deter=[];  %Give the combinations of fipi and fix that imply determinacy

%eigenvalues
eigv1=zeros(n,m);
eigv2=zeros(n,m);
eigv3=zeros(n,m);
count=zeros(n,m);

for j=1:m+1
    fix=(j-1)/decimals;
    fixs(j)=fix;
    for s=1:n+1
        fipi=(s-1)/decimals;
        fipis(s)=fipi;
        iteration=j*s

%************ Calculate jacobian and eigenvalues***************************
fx=fix+tau*(gama*fipi+fis)/(1-gama);
fpi=fipi+fis;
fs=tau*(gama*fipi+fis)/(1-gama);

B=[w+tau*si mu 0; 0 beta 0; 0 0 1];
d11=1+mu*fx; 
d12=mu*fpi;
d13=-mu*fs;
d21=-landa*k1; 
d22=1;
d23=0;
d31=1;
d32=0;
d33=0;

D=[d11 d12 d13; d21 d22 d23; d31 d32 d33];
   
f=inv(B)*D;
[eigvec,eigval]=eig(f);
egv1=eigval(1,1);
egv2=eigval(2,2);
egv3=eigval(3,3);
eigv1(s,j)=egv1;
eigv2(s,j)=egv2;
eigv3(s,j)=egv3;
counter=0;
%Let's count the number of unstable eigenvalues
    if abs(real(egv1))>1
        counter=counter+1;
    end  
    if abs(real(egv2))>1
        counter=counter+1;
    end
    if abs(real(egv3))>1
        counter=counter+1;
    end  
    count(s,j)=counter;
 %Let's analyze stability       
        if counter==2                
           deter1=[-1,-1]; 
        else
           deter1=[fipi,fix];           
        end
        deter=[deter;deter1];

    end
end

obs1=length(deter);
M1=sortrows(deter,2); 
COMB=M1(obs1,:);  %Length of indeterminacy

display('INDETERMINACY REGION FOR CPITR')
COMB
figure;
hold on

subplot(1,1,1)
plotmatrix(deter(:,1),deter(:,2))
title('Open economy with complete markets and CPITR');
xlabel('policy reaction to inflation');
ylabel('policy reaction to output-gap');

hold off
pause;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%::::::::::Forecast-based CPI inflation Taylor Rule (FB-CPITR):::::::::
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%same notation as before by including "e" at the end

%Creating output matrix
fixse=zeros(m,1); %parameter of output gap
fipise=zeros(n,1);   %parameter of inflation gap
detere=[];  %Give the combinations of fipi and fix that imply determinacy

%eigenvalues
eigve1=zeros(n,m);
eigve2=zeros(n,m);

for j=1:m+1
    fixe=(j-1)/decimals;
    fixse(j)=fixe;
    for s=1:n+1
        fipie=(s-1)/decimals;
        fipise(s)=fipie;
        iteration=j*s
              
%************ Calculate jacobian and eigenvalues***************************
fxe=fixe+tau*(gama*fipie+fis)/(1-gama);
fpie=fipie+fis;
fse=tau*(gama*fipie+fis)/(1-gama);

be11 =w-mu*fxe+tau*si;
be12 =mu*(1-fpie);
be21 = 0;
be22 = beta;
B=[be11 be12; be21 be22];

de11 =1-tau*fse;
de12 =0;
de21 = -landa*k1;
de22 = 1;
D=[de11 de12; de21 de22];
fe=inv(B)*D;
[eigvece,eigvale]=eig(fe);
egve1=eigvale(1,1);
egve2=eigvale(2,2);

eigve1(s,j)=egve1;
eigve2(s,j)=egve2;

        if abs(real(egve1))>1&&abs(real(egve2))>1               
           detere1=[-1,-1]; 
        else
           detere1=[fipie,fixe];         
        end
       detere=[detere;detere1];
 
    end
end

Me1=sortrows(detere,1);

display('INDETERMINACY REGION FOR FB-CPITR')
figure;
hold on

subplot(1,1,1)
plotmatrix(detere(:,1),detere(:,2))
title('Open economy with complete markets and FB-CPITR');
xlabel('policy reaction to inflation');
ylabel('policy reaction to output-gap');

hold off

