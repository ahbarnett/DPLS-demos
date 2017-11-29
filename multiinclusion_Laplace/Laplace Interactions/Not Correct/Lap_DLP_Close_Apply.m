function [z,zx,zy]=Lap_DLP_Close_Apply(NPt,s,t,f)

M=length(s.len);

if isfield(NPt,'s')
    z.s=zeros(length(s.x),1);
    zx.s=z.s;
    zy.s=z.s;
    ls=0;
    for k=1:M
        if ~isempty(NPt.s{k})
            ss.x=s.x(ls+1:ls+s.len(k));
            ss.t=s.t(ls+1:ls+s.len(k));
            ss.xp=s.xp(ls+1:ls+s.len(k));
            ss.sp=s.sp(ls+1:ls+s.len(k));
            ss.w=s.w(ls+1:ls+s.len(k));
            ss.nx=s.nx(ls+1:ls+s.len(k));
            ss.cw=s.cw(ls+1:ls+s.len(k));
            xx.x=s.x(NPt.s{k});
            ff=f(ls+1:ls+s.len(k));
            %in=inpoly([real(xx),imag(xx)],[real(ss.x),imag(ss.x)]);
            in=inpolygon(real(xx.x),imag(xx.x),real(ss.x),imag(ss.x));
            zz=zeros(length(xx.x),1);
            zzx=zz;
            zzy=zz;
            sin=sum(in);
            if sin>0
                [u,ux,uy]=lapDevalclose(xx.x(in),ss,ff,'i');
                zz(in)=u;
                zzx(in)=ux;
                zzy(in)=uy;
            end
            if sin<length(in)
                [u,ux,uy]=lapDevalclose(xx.x(~in),ss,ff,'e');
                zz(~in)=u;
                zzx(~in)=ux;
                zzy(~in)=uy;
            end
            DLP=LapDLPmatrix(xx,ss,0)*ff;
            xx.nx=0*xx.x+1;
            [~,DLPx]=LapDLPmatrix(xx,ss,0);
            DLPx=DLPx*ff;
            xx.nx=0*xx.x+1i;
            [~,DLPy]=LapDLPmatrix(xx,ss,0);
            DLPy=DLPy*ff;
            zz=zz-DLP;
            zzx=zzx-DLPx;
            zzy=zzy-DLPy;
            z.s(NPt.s{k})=z.s(NPt.s{k})+zz;
            zx.s(NPt.s{k})=zx.s(NPt.s{k})+zzx;
            zy.s(NPt.s{k})=zy.s(NPt.s{k})+zzy;
        end
        ls=ls+s.len(k);
    end
end

if isfield(NPt,'t')
    z.t=zeros(length(t.x),1);
    zx.t=z.t;
    zy.t=z.t;
    ls=0;
    for k=1:M
        if ~isempty(NPt.t{k})
            ss.x=s.x(ls+1:ls+s.len(k));
            ss.t=s.t(ls+1:ls+s.len(k));
            ss.xp=s.xp(ls+1:ls+s.len(k));
            ss.sp=s.sp(ls+1:ls+s.len(k));
            ss.w=s.w(ls+1:ls+s.len(k));
            ss.nx=s.nx(ls+1:ls+s.len(k));
            ss.cw=s.cw(ls+1:ls+s.len(k));
            xx.x=t.x(NPt.t{k});
            ff=f(ls+1:ls+s.len(k));
            %in=inpoly([real(xx),imag(xx)],[real(ss.x),imag(ss.x)]);
            in=inpolygon(real(xx.x),imag(xx.x),real(ss.x),imag(ss.x));
            zz=zeros(length(xx.x),1);
            zzx=zz;
            zzy=zz;
            sin=sum(in);
            if sin>0
                [u,ux,uy]=lapDevalclose(xx.x(in),ss,ff,'i');
                zz(in)=u;
                zzx(in)=ux;
                zzy(in)=uy;
            end
            if sin<length(in)
                [u,ux,uy]=lapDevalclose(xx.x(~in),ss,ff,'e');
                zz(~in)=u;
                zzx(~in)=ux;
                zzy(~in)=uy;
            end
            DLP=LapDLPmatrix(xx,ss,0)*ff;
            xx.nx=0*xx.x+1;
            [~,DLPx]=LapDLPmatrix(xx,ss,0);
            DLPx=DLPx*ff;
            xx.nx=0*xx.x+1i;
            [~,DLPy]=LapDLPmatrix(xx,ss,0);
            DLPy=DLPy*ff;
            zz=zz-DLP;
            zzx=zzx-DLPx;
            zzy=zzy-DLPy;
            z.t(NPt.t{k})=z.t(NPt.t{k})+zz;
            zx.t(NPt.t{k})=zx.t(NPt.t{k})+zzx;
            zy.t(NPt.t{k})=zy.t(NPt.t{k})+zzy;
        end
        ls=ls+s.len(k);
    end
end