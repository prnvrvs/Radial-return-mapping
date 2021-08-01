    clc;
    clear;
    % written By Pranav Kumar, IIT Madras
    % contact: prnvkmr4@gmail.com
    %====================================================
    % strain controlled combined hardening of elastic plastic system
     sigmay=400  % yield stress
     E=220e3     % Young's Modulus
    ubmax=0.010
    h=0.0001;
    ub1=0:h:ubmax;                               
    ub2=ubmax:-1*h:-1*ubmax;
    ub3=-1*ubmax:h:0;
    et=[ub1,ub2,ub3]; % strain range in one cycle
    % number of cycle 
    n=1; % number of cycle plus 1 and should be more than 1
    for k=2:1:n
        et=[et,et];
    end
    %===================================================
    S(1)=0;         % Plastic arc length
    alpha(1)=0 ;    % Back stress
    sigma(1)=0;     % initial stress
    sib(1)=0;
    K= 50e3;   % tangent modulus_kinematic
    H=10e3;     % tangent modulus_isotropic
    ep(1)=0;    % initial Plastic strain

    %====================================================
        for i=1:length(et)-1

       et_dot=et(i+1)-et(i);
       lamda_dot=0;
       T_sdot=E.*et_dot;        % incremental stress
       T_stress=sigma(i)+T_sdot;% Trial stress
       T_alpha=alpha(i);    % Trial back stress
       T_S= S(i);           % Trial plastic arc length
       f=abs(T_stress-T_alpha)-abs(sigmay+H.*lamda_dot)
       %=======================
        if T_stress-T_alpha<0
            sign_sig=-1;
        else
            sign_sig=1;
        end
        %=======================
       if f<=0
           sigma(i+1)=sigma(i)+T_sdot;
           sib(i+1)=sib(i)+T_sdot;
           alpha(i+1)=T_alpha;
           ep(i+1)=ep(i);
           S(i+1)=T_S;
           sigmay=sigmay+H.*lamda_dot
       else
           lamda_dot=f./(E+ H +K);
           ep_dot=lamda_dot.*sign_sig;
           ep(i+1)=ep(i)+ep_dot;
           alpha(i+1)=alpha(i)+K.*ep_dot;
           S(i+1)=S(i)+H.*lamda_dot;
           sib(i+1)=E.*(et(i+1)-ep(i+1));
           sigma(i+1)=T_stress-E*lamda_dot*sign_sig;
           
           sigmay=sigmay+H.*lamda_dot;
       end
           
      

        end
    %======================================================
    % animated plot
    hh=figure
       plot(nan,nan)
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
%     xlim([-ubmax-2*et_dot ubmax+2*et_dot])
%     ylim([-sigmay-100 sigmay+100])
    grid on
    grid minor
    title('Stress vs strain (combined Hardening)')
    hold on
    
    h = animatedline;
filename = 'testAnimated.gif';
    for j=1:length(et)
        addpoints(h,et(j),sigma(j));
        drawnow
            % Capture the plot as an image 
      frame = getframe(hh); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256);   pause(0.001)
       if j == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
    end
