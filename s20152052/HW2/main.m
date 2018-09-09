clear all
clc

%% N vs E_ground_err %%


%% input parameters

N=[5, 50, 500]; % point #


a=5.*1e-9; % scale [m]
m=0.19; % electron mass

param=[a; m];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(N);
    [A(1,i),A(2,i),A(3,i)]=E_ISW(N(i),param); % E_err check A(1): ground_Anl, A(2):ground_num, A(3): E_err 
end


subplot(2,1,1)
plot(N,A(1,:),'r',N,A(2,:),'b')
xlabel('discretization #')
ylabel('ground energy [eV]')
legend('Analitical','Numerical')

subplot(2,1,2)
semilogy(N,A(3,:))
xlabel('discretization #')
ylabel('ground energy error [eV]')
