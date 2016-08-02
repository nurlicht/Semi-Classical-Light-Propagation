%The source is not transform-limited.

close all
clear
clc

%Yee: Magnetic field at the edge AND after the electric field:
%Hs(m,n)==H((m+1/2)*dz,(n+1/2)*dt)

%Advance E; Advance H;

%Material properties are calculated at the cell center, and henec must
%coincide with the electric field.

%The coupling coefficient gamma is attributed to the middle of the cell.



%<Phase information>
Ph0=pi*(0);
Ph1=0;
Ph2=pi*(0);
Ph3=0;
%</Phase information>


%<Physical Constants>
h_slash=1.054571628e-34;
c=2.99792458e8;
Mu0=(4*pi)*1e-7;
Ep0=1/(Mu0*c^2);
%</Physical Constants>

%<Geometry>
N_PML=0;
N_Quantum=18000/8;
N_air=2000;
%</Geometry>

I=sqrt(-1);

%<IO Locations>
Observation_Point_1=N_PML+N_air;
Observation_Point_2=N_PML+N_air+1+fix(N_Quantum/4);
Source_Index=1;
%</IO Locations>

%<Material Properties>
N_atom=1e24;
Gamma=1e-29;
Rho30=-1;
T1=1e-10;
T2=5e-10;
%</Material Properties>

%<Source Dynamics>
f0=2e14;
Tp=100e-15;
T_delay=2*T2;
switch menu('Please select the simulation (Ref: PRA, 52, 3082)','Figure 2','Figure 3','Figure 4')
    case 1
        Pulse_Factor=1;
    case 2
        Pulse_Factor=1/2;
    case 3
        Pulse_Factor=2;
end
Emax=4.2186e9*Pulse_Factor;
Beta=1e-4;
w0=2*pi*f0;
Lambda0=c/f0;
%</Source Dynamics>

%<FDTD Parameters>
Spatial_Resolution_Factor=200;
Courant_Factor=2;
T_simulation=20*(1/f0);
dz=Lambda0/Spatial_Resolution_Factor;
dt=(dz/c)/Courant_Factor;
N_t=fix(T_simulation/dt);
N_cells=2*(N_PML+N_air)+N_Quantum;
t_center=dt*(1:N_t); %t_center(n)=n*dt, when En is calculated.
z_center=dz*(1:N_cells);  %z_center(n)=n*dz, where En is calculated.
z_center_Q=z_center(N_air+(1:N_Quantum));
Material=[0*ones(N_PML,1);1*ones(N_air,1);2*ones(N_Quantum,1);1*ones(N_air,1);0*ones(N_PML,1)];
Max_Iteration=4;
%</FDTD Parameters>

%<Source Field>
%Np=fix(Tp/dt);
%x=(2*t_center/Tp)-1;
%f_=zeros(1,N_t);
%f_(1:Np)=-4.201355*x(1:Np).*(1-x(1:Np).^2).^3;
%Ng=5*Np;
%g_=zeros(1,N_t);
%g_(1:Ng)=sin(w0*t_center(1:Ng)).*(1-x(1:Ng).^2).^4;
%if Ng<N_t
%    g_((Ng+1):N_t)=ones(1,N_t-Ng);
%end
%Delay_Index=fix(T_delay/dt);
%g_delayed=[zeros(1,Delay_Index) g_(1:(N_t-Delay_Index))];
%f_delayed=[zeros(1,Delay_Index) f_(1:(N_t-Delay_Index))];
%E_source=Emax*(f_+1*Beta*g_delayed);
%E_source=Emax*f_;
%E_source=Emax*sech(10*(t_center-Tp/2)/(Tp/2)).*sin(w0*t_center);
E_source=Emax*sech(10*(t_center-Tp/2)/(Tp/2)).*sin(w0*t_center+Ph0+Ph1*(t_center/Tp)+Ph2*(t_center/Tp).^2+Ph3*(t_center/Tp).^3);
%</Source Field>

%<Field Initialization>
E=zeros(N_cells,1);
Hs=zeros(N_cells,1);
u1=zeros(N_Quantum,1);
u2=zeros(N_Quantum,1);
u3=zeros(N_Quantum,1);
Field_1=zeros(N_t,1);
Field_2=zeros(N_t,1);
%</Field Initialization>

%<Efficient Coefficients>
Coeff1=(dt/(Mu0*dz));
Coeff2=(dt/(Ep0*dz));
Coeff3=(dt*w0/2);
%</Efficient Coefficients>

%<Efficient exponentials>
As=(N_atom*Gamma/(Ep0*T2))*exp(-(t_center+0.5*dt)/T2);
As_md=As*(1/2)*dt;
clear As;
Bs=(N_atom*Gamma*w0/Ep0)*exp(-(t_center+0.5*dt)/T2);
Bs_md=Bs*(1/2)*dt;
clear Bs;
Cs_plus=(2*Gamma/h_slash)*exp(-(t_center+0.5*dt)*(1/T1-1/T2));
Cs_plus_md=(dt*Cs_plus*(1/2)*(1/2));
clear Cs_plus;
Cs_minus=(2*Gamma/h_slash)*exp(-(t_center+0.5*dt)*(1/T2-1/T1));
Cs_minus_md=dt*Cs_minus*(1/2)*(1/2);
clear Cs_minus;
Ds=(2*Gamma/h_slash)*Rho30*exp((t_center+0.5*dt)/T2);  %New version of D
Ds_md=(1/2)*Ds*dt;
clear Ds;
%</Efficient exponentials>

%<Efficient Indices>
Index0=2:(Source_Index-1);
Index0m1=Index0-1;
Index1=1:(N_cells-1);
Index1p1=Index1+1;
Index2=2:N_air;
Index2m1=Index2-1;
Index3=N_air+(1:N_Quantum);
Index3m1=Index3-1;
Index3mOffset=Index3-min(Index3)+1;
Index4=N_air+N_Quantum+(1:N_air);
Index4m1=Index4-1;
Index5=(Source_Index+1):N_air;
Index5m1=Index5-1;
%</Efficient Indices>

U_old=zeros(N_Quantum,4);
U_new=zeros(N_Quantum,4);
Tempx=zeros(N_Quantum,1);
Tempy=zeros(N_Quantum,1);
Tempz=zeros(N_Quantum,1);

%<Time Evolution Loop>
Monitor_Index=50;
for n=1:N_t
    %<Update E>
    %<Update E in air (left)>
    E(Source_Index)=E_source(n); %HARD SOURCE
    E(Index5)=E(Index5)-Coeff2*(Hs(Index5)-Hs(Index5m1));
    %</Update E in air (left)>

    %<Update E(m) and ui(m) in Quantum medium>
    U_old=[E(Index3),u1,u2,u3];
    U_new=U_old;
    cntr=1;
    while (cntr<Max_Iteration)
%        U_=U_new;
        Tempx=(U_new(Index3mOffset,3)+U_old(Index3mOffset,3));
        Tempy=(U_new(Index3mOffset,2)+U_old(Index3mOffset,2));
        Tempz=(U_new(Index3mOffset,1)+U_old(Index3mOffset,1));
        U_new=U_old+[-Coeff2*(Hs(Index3)-Hs(Index3m1))+...
            -As_md(n)*Tempy+Bs_md(n)*Tempx,...
            +Coeff3*Tempx,...
            -Coeff3*Tempy+...
            (Cs_plus_md(n)*(U_new(Index3mOffset,4)+U_old(Index3mOffset,4))+Ds_md(n)).*Tempz,...
            -Cs_minus_md(n)*Tempz.*Tempx];
        cntr=cntr+1;
    end
    E(Index3)=U_new(:,1); 
    u1=U_new(:,2);
    u2=U_new(:,3);
    u3=U_new(:,4);
    %</Update E(m) and ui(m) in Quantum medium>

    %<Update E in air (right)>
    E(Index4)=E(Index4)-Coeff2*(Hs(Index4)-Hs(Index4m1));
    %</Update E in air (right)>
    %</Update E>

    %<Update H>
    Hs(Index1)=Hs(Index1)-Coeff1*(E(Index1p1)-E(Index1));
    Hs(N_cells)=Hs(N_cells)+Coeff1*E(N_cells);
    %</Update H>

%    Field_1(n)=exp(-t_center(n)/T2)*(u1(Observation_Point_2-N_air)+I*u2(Observation_Point_2-N_air));
    Field_1(n)=exp(-t_center(n)/T2)*u1(Observation_Point_2-N_air);
    Field_2(n)=E(Observation_Point_2);

    if n==Monitor_Index*fix(n/Monitor_Index)
        subplot(311)
        plot(z_center*1e6,E/max(abs(E)),'b',(z_center+dz/2)*1e6,Hs/max(abs(Hs)),'r',...
            1e6*(N_air+1/2)*dz*[1 1],[-1 1],'c--',1e6*(N_air+N_Quantum+1/2)*dz*[1 1],[-1 1],'c--');
        ylabel('Normalized Field');
        xlim(1e6*[z_center(1)-dz/2 z_center(N_cells)+dz/2]);
        xlabel('Position, x(\mum)');
        legend('E','H');
        title(['n=' num2str(n) ', n(Inf+NAN)=' num2str(N_Quantum*4-sum(sum(isfinite(U_new))))]);
%       title(['n=' num2str(n) ', n(Inf+NAN)=' num2str(N_Quantum*4-sum(sum(isfinite(U_new)))) ', Error=' num2str(max(max(abs(U_new-U_)./(abs(U_+1e-200)))))]);

        subplot(312)
        Temp=exp(-t_center(n)/T2);
        Rho_1=Temp*u1;
        Rho_2=Temp*u2;
        Rho_3=Rho30+u3*exp(-t_center(n)/T1);
        Q_E=Rho_1.^2+Rho_2.^2+Rho_3.^2;
        plot(z_center_Q*1e6,Rho_1,'b',z_center_Q*1e6,Rho_2,'r',z_center_Q*1e6,Rho_3,'k',z_center_Q*1e6,Q_E,'m',...
            1e6*(N_air+1/2)*dz*[1 1],[-1 1],'c--',1e6*(N_air+N_Quantum+1/2)*dz*[1 1],[-1 1],'c--');
        xlim(1e6*[z_center(1)-dz/2 z_center(N_cells)+dz/2]);
        title(['\rho1_m_a_x=' num2str(max(abs(Rho_1))) ', \rho2_m_a_x=' num2str(max(abs(Rho_2))) ...
            ', \rho3_m_a_x=' num2str(max(abs(Rho_3))) ', \Sigma\rho^2=' num2str(min(Q_E)) ':' num2str(max(Q_E))]);
        xlabel('Position, x(\mum)');
        legend('\rho1','\rho2','\rho3','\Sigma\rho^2');

        subplot(313)
        Temp=dt*(1:n)*1e15;
        plot(Temp,Field_1(1:n)/max(abs(Field_1(1:n))),'b',Temp,Field_2(1:n)/max(abs(Field_2(1:n))),'r',n*dt*1e15*[1 1],[0 1],'--k');
        legend('\rho_1(t)','E(t)');
        xlabel('Time, t(fs)');
        ylabel('Observed Field');
%        xlim(1e15*[0 N_cells*dz/c])
%        xlim(1e15*dt*[1 N_t])
        xlim([00 200])

        pause(0.001);
    end
end
%</Time Evolution Loop>

f_max_THz=1e3;
N_new=fix(N_t/1);
NN=50;

Index_trun=1:N_new;
%Window=hanning(N_new);
Window=1;
Rho1_f=fft(real(Field_1(Index_trun)).*Window);
E_f=fft(real(Field_2(Index_trun)).*Window);
Refractive_Index=sqrt(1-N_atom*Gamma*Rho1_f./(Ep0*E_f));
Xi=-N_atom*Gamma*Rho1_f./(Ep0*E_f);
FFF=1e-12*(Index_trun-1)/(N_new*dt);
cntr=1;
while FFF(cntr)<f_max_THz
    cntr=cntr+1;
end
Nf_disp=cntr;
FFF=FFF(1:Nf_disp);
Xi=Xi(1:Nf_disp);
Refractive_Index=Refractive_Index(1:Nf_disp);
%plot(FFF,10*filtfilt(ones(NN,1)/NN,1,imag(Refractive_Index)),'b-',FFF,filtfilt(ones(NN,1)/NN,1,real(Refractive_Index)),'r',FFF,E_f/max(abs(E_f(Index_trun))),'k');
subplot(311)
title(['\phi_0=' num2str(Ph0) ', \phi_1=' num2str(Ph1) ', \phi_2=' num2str(Ph2) ', \phi_3=' num2str(Ph3) ', Pulse: ' num2str(Pulse_Factor*2) '\pi'])

figure
subplot(221)
Temp=filtfilt(ones(NN,1)/NN,1,imag(Xi));
plot(FFF,Temp,'b-');
ylabel('Imag(\chi)');
xlabel('Frequency, f (THz)');
xlim([0 f_max_THz]);
ylim([min(Temp) max(Temp)]);
title(['\phi_0=' num2str(Ph0) ', \phi_1=' num2str(Ph1) ', \phi_2=' num2str(Ph2) ', \phi_3=' num2str(Ph3)])
subplot(223)
Temp=filtfilt(ones(NN,1)/NN,1,real(Xi));
plot(FFF,Temp,'r');
xlabel('Frequency, f (THz)');
ylabel('Real(\chi)');
ylim([min(Temp) max(Temp)]);
xlim([0 f_max_THz]);

subplot(222)
Temp=filtfilt(ones(NN,1)/NN,1,imag(Refractive_Index));
plot(FFF,Temp,'b-');
ylabel('Absorption');
xlabel('Frequency, f (THz)');
xlim([0 f_max_THz]);
ylim([min(Temp) max(Temp)]);
title(['Pulse: ' num2str(Pulse_Factor*2) '\pi'])
subplot(224)
Temp=filtfilt(ones(NN,1)/NN,1,real(Refractive_Index));
plot(FFF,Temp,'r');
xlabel('Frequency, f (THz)');
ylabel('Refractive Index');
ylim([min(Temp) max(Temp)]);
xlim([0 f_max_THz]);
