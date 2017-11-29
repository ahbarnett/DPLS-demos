function u=SLP_Far_Apply(s,if_s,t,if_t,f,iprec_S)

if nargin==5
    iprec_S=4;
end

if isfloat(iprec_S)
    tau = real(f);
    [I1x1, I3x1, I3x2] = SLP_Lap_FMM(s,if_s,t,if_t,tau,iprec_S);
    tau = imag(f);
    [I1x2, I4x1, I4x2] = SLP_Lap_FMM(s,if_s,t,if_t,tau,iprec_S);
    tau = real(s.x.*conj(f));
    [~, I2x1, I2x2] = SLP_Lap_FMM(s,if_s,t,if_t,tau,iprec_S);
    
    if if_s==1
        I1.s = (I1x1.s+1i*I1x2.s)/2;
        I2.s = (I2x1.s+1i*I2x2.s)/2;
        I3.s = real(s.x).*(I3x1.s+1i*I3x2.s)/2;
        I4.s = imag(s.x).*(I4x1.s+1i*I4x2.s)/2;
        u.s = I1.s+I2.s-I3.s-I4.s;
    end
    
    if if_t==1
        I1.t = (I1x1.t+1i*I1x2.t)/2;
        I2.t = (I2x1.t+1i*I2x2.t)/2;
        I3.t = real(t.x).*(I3x1.t+1i*I3x2.t)/2;
        I4.t = imag(t.x).*(I4x1.t+1i*I4x2.t)/2;
        u.t = I1.t+I2.t-I3.t-I4.t;
    end
else
    if if_s==1
        u.s=iprec_S.s*[real(f);imag(f)];
        u.s=u.s(1:end/2)+1i*u.s(end/2+1:end);
    end
    if if_t==1
        u.t=iprec_S.t*[real(f);imag(f)];
        u.t=u.t(1:end/2)+1i*u.t(end/2+1:end);
    end
end
end