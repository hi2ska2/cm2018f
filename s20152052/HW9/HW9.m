clear all
clc

eps_si = 11.7; eps_ox = 3.9; % Relative permittivity 
ny=9;
nz=5;
Nnode=ny*nz;

y=[1:ny];
z=[1:nz];
A=zeros(Nnode,Nnode);

% %% case1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % empty circle
% for jj=2:nz-1 % z-axis 
%     for ii=2:ny-1 % y-axis
% %         ind=ii+ny*(jj-1); %index
%         A(ii+ny*(jj-1),ii+1+ny*(jj-1))=1; 
%         A(ii+ny*(jj-1),ii+ny*(jj))=1; 
%         A(ii+ny*(jj-1),ii+ny*(jj-1))=-4; 
%         A(ii+ny*(jj-1),ii-1+ny*(jj-1))=1;
%         A(ii+ny*(jj-1),ii+ny*(jj-2))=1;
%     end
% end
% 
% % black circle
% for jj=2:4
%     for ii=1
%         A(ii+ny*(jj-1),ii+1+ny*(jj-1))=1; 
%         A(ii+ny*(jj-1),ii+ny*(jj))=0.5; 
%         A(ii+ny*(jj-1),ii+ny*(jj-1))=-2; 
%         A(ii+ny*(jj-1),ii+ny*(jj-2))=0.5; 
%     end
%     for ii=9
%         A(ii+ny*(jj-1),ii+ny*(jj))=0.5; 
%         A(ii+ny*(jj-1),ii+ny*(jj-1))=-2; 
%         A(ii+ny*(jj-1),ii+ny*(jj-2))=0.5; 
%         A(ii+ny*(jj-1),ii-1+ny*(jj-1))=1; 
%     end
% end
% for jj=5
%     for ii=3:7
%         A(ii+ny*(jj-1),ii+1+ny*(jj-1))=0.5; 
%         A(ii+ny*(jj-1),ii+ny*(jj-1))=-2; 
%         A(ii+ny*(jj-1),ii+ny*(jj-2))=1; 
%         A(ii+ny*(jj-1),ii-1+ny*(jj-1))=0.5; 
%     end
% end
% 
% % blue circle
% for jj=1
%     for ii=1:ny
%         A(ii+ny*(jj-1),ii+ny*(jj-1))=1;
%     end
% end
% 
% for jj=5
%     for ii=1:2
%         A(ii+ny*(jj-1),ii+ny*(jj-1))=1;
%     end
% end
% 
% % red circle
% for jj=5
%     for ii=8:9
%         A(ii+ny*(jj-1),ii+ny*(jj-1))=1;
%     end
% end
% b=zeros(Nnode,1);
% b(Nnode,1)=1;
% b(Nnode-1,1)=1;
% 
% 
% psi=A\b;
% 
% psi_2d=reshape(psi,ny,nz);
% 
% surf(psi_2d)
% xlabel('z-axis')
% ylabel('y-axis')
% zlabel('psi')
% 

% %% case2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % empty circle
% for jj=2:nz-1 % z-axis 
%     for ii=2:ny-1 % y-axis
% %         ind=ii+ny*(jj-1); %index
%         A(ii+ny*(jj-1),ii+1+ny*(jj-1))=1; 
%         A(ii+ny*(jj-1),ii+ny*(jj))=1; 
%         A(ii+ny*(jj-1),ii+ny*(jj-1))=-4; 
%         A(ii+ny*(jj-1),ii-1+ny*(jj-1))=1;
%         A(ii+ny*(jj-1),ii+ny*(jj-2))=1;
%     end
% end
% 
% % black circle
% for jj=2:4
%     for ii=9
%         A(ii+ny*(jj-1),ii+ny*(jj))=0.5; 
%         A(ii+ny*(jj-1),ii+ny*(jj-1))=-2; 
%         A(ii+ny*(jj-1),ii+ny*(jj-2))=0.5; 
%         A(ii+ny*(jj-1),ii-1+ny*(jj-1))=1; 
%     end
% end
% 
% % blue circle
% for jj=1
%     for ii=1:ny
%         A(ii+ny*(jj-1),ii+ny*(jj-1))=1;
%     end
% end
% 
% % red circle
% for jj=5
%     for ii=1:9
%         A(ii+ny*(jj-1),ii+ny*(jj-1))=1;
%     end
% end
% for jj=2:4
%     for ii=1
%         A(ii+ny*(jj-1),ii+ny*(jj-1))=1;
%     end
% end
% 
% b=zeros(Nnode,1);
% b(37:45,1)=1;
% b(10,1)=1;
% b(19,1)=1;
% b(28,1)=1;
% 
% 
% psi=A\b;
% 
% psi_2d=reshape(psi,ny,nz);
% 
% surf(psi_2d)
% xlabel('z-axis')
% ylabel('y-axis')
% zlabel('psi')




% %% case3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % empty circle
% for jj=2:nz-1 % z-axis 
%     for ii=2:ny-1 % y-axis
% %         ind=ii+ny*(jj-1); %index
%         A(ii+ny*(jj-1),ii+1+ny*(jj-1))=1; 
%         A(ii+ny*(jj-1),ii+ny*(jj))=1; 
%         A(ii+ny*(jj-1),ii+ny*(jj-1))=-4; 
%         A(ii+ny*(jj-1),ii-1+ny*(jj-1))=1;
%         A(ii+ny*(jj-1),ii+ny*(jj-2))=1;
%     end
% end
% 
% % black circle
% for jj=2:4
%     for ii=1
%         A(ii+ny*(jj-1),ii+1+ny*(jj-1))=1; 
%         A(ii+ny*(jj-1),ii+ny*(jj))=0.5; 
%         A(ii+ny*(jj-1),ii+ny*(jj-1))=-2; 
%         A(ii+ny*(jj-1),ii+ny*(jj-2))=0.5; 
%     end
%     for ii=9
%         A(ii+ny*(jj-1),ii+ny*(jj))=0.5; 
%         A(ii+ny*(jj-1),ii+ny*(jj-1))=-2; 
%         A(ii+ny*(jj-1),ii+ny*(jj-2))=0.5; 
%         A(ii+ny*(jj-1),ii-1+ny*(jj-1))=1; 
%     end
% end
% for jj=5
%     for ii=3:7
%         A(ii+ny*(jj-1),ii+1+ny*(jj-1))=0.5; 
%         A(ii+ny*(jj-1),ii+ny*(jj-1))=-2; 
%         A(ii+ny*(jj-1),ii+ny*(jj-2))=1; 
%         A(ii+ny*(jj-1),ii-1+ny*(jj-1))=0.5; 
%     end
% end
% 
% % blue circle
% for jj=5
%     for ii=1:2
%         A(ii+ny*(jj-1),ii+ny*(jj-1))=1;
%     end
% end
% 
% % red circle
% for jj=5
%     for ii=8:9
%         A(ii+ny*(jj-1),ii+ny*(jj-1))=1;
%     end
% end
% for jj=1
%     for ii=1:9
%         A(ii+ny*(jj-1),ii+ny*(jj-1))=1;
%     end
% end
% 
% b=zeros(Nnode,1);
% b(1:9,1)=1;
% b(Nnode,1)=1;
% b(Nnode-1,1)=1;
% 
% 
% psi=A\b;
% 
% psi_2d=reshape(psi,ny,nz);
% 
% surf(psi_2d)
% xlabel('z-axis')
% ylabel('y-axis')
% zlabel('psi')



%% case4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% empty circle
for jj=2:nz-1 % z-axis 
    for ii=2:ny-1 % y-axis
%         ind=ii+ny*(jj-1); %index
        A(ii+ny*(jj-1),ii+1+ny*(jj-1))=1; 
        A(ii+ny*(jj-1),ii+ny*(jj))=1; 
        A(ii+ny*(jj-1),ii+ny*(jj-1))=-4; 
        A(ii+ny*(jj-1),ii-1+ny*(jj-1))=1;
        A(ii+ny*(jj-1),ii+ny*(jj-2))=1;
    end
end

% black circle
for jj=2:4
    for ii=9
        A(ii+ny*(jj-1),ii+ny*(jj))=0.5; 
        A(ii+ny*(jj-1),ii+ny*(jj-1))=-2; 
        A(ii+ny*(jj-1),ii+ny*(jj-2))=0.5; 
        A(ii+ny*(jj-1),ii-1+ny*(jj-1))=1; 
    end
end

% red circle
for jj=5
    for ii=1:9
        A(ii+ny*(jj-1),ii+ny*(jj-1))=1;
    end
end

for jj=1
    for ii=1:9
        A(ii+ny*(jj-1),ii+ny*(jj-1))=1;
    end
end
for jj=2:4
    for ii=1
        A(ii+ny*(jj-1),ii+ny*(jj-1))=1;
    end
end

b=zeros(Nnode,1);
b(1:9,1)=1;
b(37:45,1)=1;
b(10,1)=1;
b(19,1)=1;
b(28,1)=1;

psi=A\b;

psi_2d=reshape(psi,ny,nz);

surf(psi_2d)
xlabel('z-axis')
ylabel('y-axis')
zlabel('psi')

