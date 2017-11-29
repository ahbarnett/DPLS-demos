function u=DLP_Far_Apply(s,if_s,t,if_t,f,iprec_D)

if nargin==5
    iprec_D=4;
end

if isfloat(iprec_D)
    tauf = f./s.nx.*real(s.nx);  % feed complex tau to Laplace close eval
    I1x1 = DLP_Lap_FMM(s,if_s,t,if_t,tauf,iprec_D);
    tauf = f./s.nx.*imag(s.nx);
    I1x2 = DLP_Lap_FMM(s,if_s,t,if_t,tauf,iprec_D);
    
    % find I_2
    tau = real(s.x.*conj(f));
    [~, I2x1, I2x2] = DLP_Lap_FMM(s,if_s,t,if_t,tau,iprec_D);
    
    % find I_3
    tau = real(f);
    [~, I3x1, I3x2] = DLP_Lap_FMM(s,if_s,t,if_t,tau,iprec_D);
    
    % find I_4
    tau = imag(f);
    [~, I4x1, I4x2] = DLP_Lap_FMM(s,if_s,t,if_t,tau,iprec_D);
    
    if if_s==1
        I1.s = I1x1.s+1i*I1x2.s;
        I2.s = I2x1.s+1i*I2x2.s;
        I3.s = real(s.x).*(I3x1.s+1i*I3x2.s);
        I4.s = imag(s.x).*(I4x1.s+1i*I4x2.s);
        u.s = I1.s+I2.s-I3.s-I4.s;
    end
    if if_t==1
        I1.t = I1x1.t+1i*I1x2.t;
        I2.t = I2x1.t+1i*I2x2.t;
        I3.t = real(t.x).*(I3x1.t+1i*I3x2.t);
        I4.t = imag(t.x).*(I4x1.t+1i*I4x2.t);
        u.t = I1.t+I2.t-I3.t-I4.t;
    end

    if if_s==1
        c = -s.cur/2/pi;           % diagonal limit of Laplace DLP
        sx = 1i*s.nx; s1=real(sx); s2=imag(sx);     % tangent vectors on src curve
        u.s=u.s+(c.*s1.^2.*s.w.*real(f)+c.*s1.*s2.*s.w.*imag(f))+1i*(c.*s1.*s2.*s.w.*real(f)+c.*s2.^2.*s.w.*imag(f));
    end
else
    if if_s==1
        u.s=iprec_D.s*[real(f);imag(f)];
        u.s=u.s(1:end/2)+1i*u.s(end/2+1:end);
    end
    if if_t==1
        u.t=iprec_D.t*[real(f);imag(f)];
        u.t=u.t(1:end/2)+1i*u.t(end/2+1:end);
    end
end
end