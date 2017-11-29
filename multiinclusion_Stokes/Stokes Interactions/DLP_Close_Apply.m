function z=DLP_Close_Apply(NPt,s,t,f,CDLP)

M=length(s.len);

if isfield(NPt,'s')
    z.s=zeros(length(s.x),1);
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
            sin=sum(in);
            if sin>0
                zz(in)=StokesDcloseeval(xx.x(in),ss,ff,'i');
            end
            if sin<length(in)
                if nargin==5
                    zz(~in)=StokesDcloseeval(xx.x(~in),ss,ff,'e',CDLP{k});
                else
                    zz(~in)=StokesDcloseeval(xx.x(~in),ss,ff,'e');
                end
            end
            DLP=DLPmatrix(xx,ss,1)*[real(ff);imag(ff)];
            zz=zz-DLP(1:end/2)-1i*DLP(end/2+1:end);
            z.s(NPt.s{k})=z.s(NPt.s{k})+zz;
        end
        ls=ls+s.len(k);
    end
end

if isfield(NPt,'t')
    z.t=zeros(length(t.x),1);
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
            sin=sum(in);
            if sin>0
                zz(in)=StokesDcloseeval(xx.x(in),ss,ff,'i');
            end
            if sin<length(in)
                if nargin==5
                    zz(~in)=StokesDcloseeval(xx.x(~in),ss,ff,'e',CDLP{k});
                else
                    zz(~in)=StokesDcloseeval(xx.x(~in),ss,ff,'e');
                end
            end
            DLP=DLPmatrix(xx,ss,1)*[real(ff);imag(ff)];
            zz=zz-DLP(1:end/2)-1i*DLP(end/2+1:end);
            z.t(NPt.t{k})=z.t(NPt.t{k})+zz;
        end
        ls=ls+s.len(k);
    end
end