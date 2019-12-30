% leap frog scheme for du/dt+du/dx =0 
% You can see diffusion error through upwind scheme

clear all;
close all;
clc;

n = 800; % number of grid points
% m = 100; % number of time steps
x = linspace(-pi,pi,n);
dx = x(2)-x(1);
CFL = 0.9;
dt = dx*CFL;
m = round(2*pi/dt);% number of time steps for 2pi 


%%  Initial condition u = sin(x)

x = linspace(-pi,pi,n);
x = [-pi-pi/n,x,pi+pi/n];
u_inital = sin(x);
u_total = zeros(m,n+2);%Two extra points on left and right
u_total(1,:) = u_inital; %Appending IC t=0

%Appending IC at t=dt UPWIND
for j = 2:n+2 % grid points
    u_total(2,j) = u_total(1,j)*(1-CFL)+u_total(1,j-1)*CFL;
end
u_total(2,1) = u_total(2,n);% Upwind cannot calculate the value at initial grid point
                            % Symmetry BC is applied

%Leap frog scheme

for i = 3:m %time
    for j = 2:n+1 % grid points
        u_total(i,j) = u_total(i-2,j)-CFL*u_total(i-1,j+1)+CFL*u_total(i-1,j-1);
    end
    u_total(i,n+2) = u_total(i,2);% Leapfrog cannot calculate the value at initial and final grid points
    u_total(i,1) = u_total(i,n); % since it is central diff scheme, Symmetry BC is applied
    disp(i)
end



%%  Initial condition u(x,0) = u(x+2kpi,0) Square wave
x = linspace(-pi,pi,n);
u_inital = zeros(1,n);
u_inital(abs(x) <= pi/2) =1; %Square wave
%Appending IC at t=dt
u_total = zeros(m,n);
u_total(1,:) = u_inital; %Appending IC

for j = 2:n % grid points
    u_total(2,j) = u_total(1,j)*(1-CFL)+u_total(1,j-1)*CFL;
end
u_total(2,1) = u_total(2,n);



%Leap frog scheme

for i = 3:m %time
    for j = 2:n-1 % grid points
        u_total(i,j) = u_total(i-2,j)-CFL*u_total(i-1,j+1)+CFL*u_total(i-1,j-1);
    end
    u_total(i,n) = u_total(i-1,n)*(1-CFL)+u_total(i-1,n-1)*CFL; %upwind for last grid point
    u_total(i,1) = u_total(i,n);
    disp(i)
end

%% Video maker


writerObj = VideoWriter(['Leapfrog_Square_pi|',num2str(n),'_CFL=',num2str(CFL),'-2'],'MPEG-4');
writerObj.FrameRate = 100;

open(writerObj);

for i=1:m    
    plot(x,u_total(i,:),'r-'); 
    %change n manually in the title everytime!
    title('Solution using leapfrog scheme for square IC; $$ \Delta x = \frac{\pi}{800}$$','interpreter','latex','FontSize',15)
    xlabel('{\it{x}} position','FontSize',15)
    ylabel('{\it{u}}','FontSize',15)
    txt = ['CFL = ',num2str(CFL)];
    text(-0,1.25,txt,'FontSize',15,'EdgeColor','black')
    axis([-pi pi -0.5 1.5])

    frame =  getframe(gcf);
    writeVideo(writerObj, frame);
end

close(writerObj);

%% Actual solution

%AFter 2pi timeperiod
act_sol =  sin(x);
act_sol2 = zeros(1,n);
act_sol2(abs(x) <= pi/2) =1; %Square wave
plot(x,act_sol2,x,u_total(end,:))

title('Solution using leapfrog scheme for square IC; $$ \Delta x = \frac{\pi}{800}$$','interpreter','latex','FontSize',15)
xlabel('{\it{x}} position','FontSize',15)
ylabel('{\it{u}}','FontSize',15)
txt = ['CFL = ',num2str(CFL)];
text(-pi/8,1.3,txt,'FontSize',15,'EdgeColor','black')
axis([-pi pi -0.5 1.5])
legend('Analytical solution','Leapfrog scheme solution')