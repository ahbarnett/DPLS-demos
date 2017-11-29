function z=Lap_SLP_Close_Apply(NPt,s,t,f)

M=length(s.len);

if isfield(NPt,'s')
    z.s=zeros(length(s.x),1);
    z.sn=z.s;
    %zx.s=z.s;
    %zy.s=z.s;
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
            xx.nx=s.nx(NPt.s{k});
            ff=f(ls+1:ls+s.len(k));
            %in=inpoly([real(xx),imag(xx)],[real(ss.x),imag(ss.x)]);
            in=inpolygon(real(xx.x),imag(xx.x),real(ss.x),imag(ss.x));
            zz=zeros(length(xx.x),1);
            zzn=zz;
%             zzx=zz;
%             zzy=zz;
            sin=sum(in);
            if sin>0
                [u,ux,uy]=lapSevalclose(xx.x(in),ss,ff,'i');
                zz(in)=u;
                zzn(in)=real(xx.nx(in)).*ux+imag(xx.nx(in)).*uy; %<----------Added Normals need to finish rest!!!
            end
            if sin<length(in)
                [u,ux,uy]=lapSevalclose(xx.x(~in),ss,ff,'e');
                zz(~in)=u;
                zzn(~in)=real(xx.nx(~in)).*ux+imag(xx.nx(~in)).*uy;
            end
            [SLP,SLPn]=LapSLPmatrix(xx,ss,0);
%             xx.nx=0*xx.x+1;
%             [~,SLPx]=LapSLPmatrix(xx,ss,0);
%             SLPx=SLPx*ff;
%             xx.nx=0*xx.x+1i;
%             [~,SLPy]=LapSLPmatrix(xx,ss,0);
%             SLPy=SLPy*ff;
            zz=zz-SLP*ff;
            zzn=zzn-SLPn*ff;
%             zzx=zzx-SLPx;
%             zzy=zzy-SLPy;
            z.s(NPt.s{k})=z.s(NPt.s{k})+zz;
            z.sn(NPt.s{k})=z.sn(NPt.s{k})+zzn;
%             zx.s(NPt.s{k})=zx.s(NPt.s{k})+zzx;
%             zy.s(NPt.s{k})=zy.s(NPt.s{k})+zzy;
        end
        ls=ls+s.len(k);
    end
end

if isfield(NPt,'t')
    z.t=zeros(length(t.x),1);
    z.tx=z.t;
    z.ty=z.t;
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
                [u,ux,uy]=lapSevalclose(xx.x(in),ss,ff,'i');
                zz(in)=u;
                zzx(in)=ux;
                zzy(in)=uy;
            end
            if sin<length(in)
                [u,ux,uy]=lapSevalclose(xx.x(~in),ss,ff,'e');
                zz(~in)=u;
                zzx(~in)=ux;
                zzy(~in)=uy;
            end
            SLP=LapSLPmatrix(xx,ss,0)*ff;
            xx.nx=0*xx.x+1;
            [~,SLPx]=LapSLPmatrix(xx,ss,0);
            SLPx=SLPx*ff;
            xx.nx=0*xx.x+1i;
            [~,SLPy]=LapSLPmatrix(xx,ss,0);
            SLPy=SLPy*ff;
            zz=zz-SLP;
            zzx=zzx-SLPx;
            zzy=zzy-SLPy;
            z.t(NPt.t{k})=z.t(NPt.t{k})+zz;
            z.tx(NPt.t{k})=z.tx(NPt.t{k})+zzx;
            z.ty(NPt.t{k})=z.ty(NPt.t{k})+zzy;
        end
        ls=ls+s.len(k);
    end
end