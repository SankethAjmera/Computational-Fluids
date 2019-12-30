% Upwind scheme for du/dt+du/dx =0 
% You can see diffusion error through upwind scheme

clear;
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

u_inital = sin(x);
u_total = zeros(m,n); %Rows are marching in time; columns are grid points in space
u_total(1,:) = u_inital; %Appending IC

for i = 2:m %time steps
    for j = 2:n % grid points
        u_total(i,j) = u_total(i-1,j)*(1-CFL)+u_total(i-1,j-1)*CFL;
    end
    u_total(i,1) = u_total(i,n);
end


%%  Initial condition u(x,0) = u(x+2kpi,0) Square wave
x = linspace(-pi,pi,n);

u_inital = zeros(1,n);
u_inital(abs(x) <= pi/2) =1; %Square wave

u_total = zeros(m,n);
u_total(1,:) = u_inital; %Appending IC
% x = [x(1)-dx,x];
for i = 2:m %time steps
    for j = 2:n % grid points
        u_total(i,j) = u_total(i-1,j)*(1-CFL)+u_total(i-1,j-1)*CFL;
    end
      u_total(i,1) = u_total(i,n);
end


%% Video maker


writerObj = VideoWriter(['Upwind_Square_pi|',num2str(n),'_CFL=',num2str(CFL),'-2'],'MPEG-4');
writerObj.FrameRate = 100;

open(writerObj);

for i=1:m    
    plot(x,u_total(i,:),'r-'); 
    
    title('Solution using upwind scheme for square IC; $$ \Delta x = \frac{\pi}{800}$$','interpreter','latex','FontSize',15)
    xlabel('{\it{x}} position','FontSize',15)
    ylabel('{\it{u}}','FontSize',15)
    txt = ['CFL = ',num2str(CFL)];
    text(-pi/8,1.1,txt,'FontSize',15,'EdgeColor','black')
    axis([-pi pi -.05 1.2])

    frame =  getframe(gcf);
    writeVideo(writerObj, frame);
end

close(writerObj);

x%% Actual solution

%AFter 2pi timeperiod
act_sol =  sin(x);
act_sol2 = zeros(1,n);
act_sol2(abs(x) <= pi/2) =1; %Square wave
plot(x,act_sol2,x,u_total(end,:))

title('Solution using upwind scheme for square IC; $$ \Delta x = \frac{\pi}{800}$$','interpreter','latex','FontSize',15)
xlabel('{\it{x}} position','FontSize',15)
ylabel('{\it{u}}','FontSize',15)
txt = ['CFL = ',num2str(CFL)];
text(-pi/8,0.7,txt,'FontSize',15,'EdgeColor','black')
axis([-pi pi -.1 1.1])
legend('Analytical solution','Upwind scheme solution')

