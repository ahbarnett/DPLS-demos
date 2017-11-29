function t=IslandIntersect(t,s,n)

if ~isempty(s)&&~isempty(t)
    sM=round(length(s)/n);
    tM=round(length(t)/n);
    
    edge=reshape((1:sM*n)',n,sM);
    edge=[(1:sM*n)',reshape(edge([2:end,1],:),sM*n,1)];
    in=inpoly(t,s,edge);
    
    in=reshape(in,n,tM);
    in=ones(n,1)*max(in);
    in=reshape(in,n*tM,1);
    t=t(~in,:);
    
    
    % Top Source, Bottom Target
    tM=round(length(t)/n);
    
    ss=reshape(s(:,2),n,sM);
    ss=max(ss);
    ss_max=max(ss);
    ss=reshape(ones(n,1)*(ss>2*pi),sM*n,1);
    ss=logical(ss);
    
    tt=reshape(t(:,2),n,tM);
    tt=min(tt);
    tt=reshape(ones(n,1)*(tt<=ss_max-2*pi),tM*n,1);
    tt=logical(tt);
    
    
    sN=round(sum(ss)/n); tN=round(sum(tt)/n);
    if sN>0&&tN>0
    edge=reshape((1:sN*n)',n,sN);
    edge=[(1:sN*n)',reshape(edge([2:end,1],:),sN*n,1)];
    in_t=inpoly(t(tt,:)+ones(tN*n,1)*[0 2*pi],s(ss,:),edge);
    in_t=reshape(in_t,n,tN);
    in_t=ones(n,1)*max(in_t);
    in_t=reshape(in_t,n*tN,1);
    in=zeros(n*tM,1);
    in(tt)=in_t;
    t=t(~in,:);
    end
    
    % Bottom Source, Top Target
    tM=round(length(t)/n);
    
    ss=reshape(s(:,2),n,sM);
    ss=min(ss);
    ss_min=min(ss);
    ss=reshape(ones(n,1)*(ss<=0),sM*n,1);
    ss=logical(ss);
    
    tt=reshape(t(:,2),n,tM);
    tt=max(tt);
    tt=reshape(ones(n,1)*(tt>ss_min+2*pi),tM*n,1);
    tt=logical(tt);
    
    
    sN=round(sum(ss)/n); tN=round(sum(tt)/n);
    if sN>0&&tN>0
    edge=reshape((1:sN*n)',n,sN);
    edge=[(1:sN*n)',reshape(edge([2:end,1],:),sN*n,1)];
    in_t=inpoly(t(tt,:)+ones(tN*n,1)*[0 -2*pi],s(ss,:),edge);
    in_t=reshape(in_t,n,tN);
    in_t=ones(n,1)*max(in_t);
    in_t=reshape(in_t,n*tN,1);
    in=zeros(n*tM,1);
    in(tt)=in_t;
    t=t(~in,:);
    end
    
    % Right Source, Left Target
    tM=round(length(t)/n);
    
    ss=reshape(s(:,1),n,sM);
    ss=max(ss);
    ss_max=max(ss);
    ss=reshape(ones(n,1)*(ss>2*pi),sM*n,1);
    ss=logical(ss);
    
    tt=reshape(t(:,1),n,tM);
    tt=min(tt);
    tt=reshape(ones(n,1)*(tt<=ss_max-2*pi),tM*n,1);
    tt=logical(tt);
    
    
    sN=round(sum(ss)/n); tN=round(sum(tt)/n);
    if sN>0&&tN>0
    edge=reshape((1:sN*n)',n,sN);
    edge=[(1:sN*n)',reshape(edge([2:end,1],:),sN*n,1)];
    in_t=inpoly(t(tt,:)+ones(tN*n,1)*[2*pi 0],s(ss,:),edge);
    in_t=reshape(in_t,n,tN);
    in_t=ones(n,1)*max(in_t);
    in_t=reshape(in_t,n*tN,1);
    in=zeros(n*tM,1);
    in(tt)=in_t;
    t=t(~in,:);
    end
    
    % Left Source, Right Target
    tM=round(length(t)/n);
    
    ss=reshape(s(:,1),n,sM);
    ss=min(ss);
    ss_min=min(ss);
    ss=reshape(ones(n,1)*(ss<=0),sM*n,1);
    ss=logical(ss);
    
    tt=reshape(t(:,1),n,tM);
    tt=max(tt);
    tt=reshape(ones(n,1)*(tt>ss_min+2*pi),tM*n,1);
    tt=logical(tt);
    
    
    sN=round(sum(ss)/n); tN=round(sum(tt)/n);
    if sN>0&&tN>0
    edge=reshape((1:sN*n)',n,sN);
    edge=[(1:sN*n)',reshape(edge([2:end,1],:),sN*n,1)];
    in_t=inpoly(t(tt,:)+ones(tN*n,1)*[-2*pi 0],s(ss,:),edge);
    in_t=reshape(in_t,n,tN);
    in_t=ones(n,1)*max(in_t);
    in_t=reshape(in_t,n*tN,1);
    in=zeros(n*tM,1);
    in(tt)=in_t;
    t=t(~in,:);
    end
    
    % Top Right Source, Bottom Left Target
    tM=round(length(t)/n);
    
    ssx=reshape(s(:,1),n,sM);
    ssx=max(ssx);
    ssx_max=max(ssx);
    ssy=reshape(s(:,2),n,sM);
    ssy=max(ssy);
    ssy_max=max(ssy);
    ss=reshape(ones(n,1)*(ssx>2*pi&ssy>2*pi),sM*n,1);
    ss=logical(ss);
    
    ttx=reshape(t(:,1),n,tM);
    ttx=min(ttx);
    tty=reshape(t(:,2),n,tM);
    tty=min(tty);
    tt=reshape(ones(n,1)*(ttx<=ssx_max-2*pi&tty<=ssy_max-2*pi),tM*n,1);
    tt=logical(tt);
    
    
    sN=round(sum(ss)/n); tN=round(sum(tt)/n);
    if sN>0&&tN>0
    edge=reshape((1:sN*n)',n,sN);
    edge=[(1:sN*n)',reshape(edge([2:end,1],:),sN*n,1)];
    in_t=inpoly(t(tt,:)+ones(tN*n,1)*[2*pi 2*pi],s(ss,:),edge);
    in_t=reshape(in_t,n,tN);
    in_t=ones(n,1)*max(in_t);
    in_t=reshape(in_t,n*tN,1);
    in=zeros(n*tM,1);
    in(tt)=in_t;
    t=t(~in,:);
    end

    % Bottom Left Source, Top Right Target
    tM=round(length(t)/n);
    
    ssx=reshape(s(:,1),n,sM);
    ssx=min(ssx);
    ssx_min=min(ssx);
    ssy=reshape(s(:,2),n,sM);
    ssy=min(ssy);
    ssy_min=min(ssy);
    ss=reshape(ones(n,1)*(ssx<=0&ssy<=0),sM*n,1);
    ss=logical(ss);
    
    ttx=reshape(t(:,1),n,tM);
    ttx=max(ttx);
    tty=reshape(t(:,2),n,tM);
    tty=max(tty);
    tt=reshape(ones(n,1)*(ttx>ssx_min+2*pi&tty>ssy_min+2*pi),tM*n,1);
    tt=logical(tt);
    
    
    sN=round(sum(ss)/n); tN=round(sum(tt)/n);
    if sN>0&&tN>0
    edge=reshape((1:sN*n)',n,sN);
    edge=[(1:sN*n)',reshape(edge([2:end,1],:),sN*n,1)];
    in_t=inpoly(t(tt,:)+ones(tN*n,1)*[-2*pi -2*pi],s(ss,:),edge);
    in_t=reshape(in_t,n,tN);
    in_t=ones(n,1)*max(in_t);
    in_t=reshape(in_t,n*tN,1);
    in=zeros(n*tM,1);
    in(tt)=in_t;
    t=t(~in,:);
    end
    
        % Top Left Source, Bottom Right Target
    tM=round(length(t)/n);
    
    ssx=reshape(s(:,1),n,sM);
    ssx=min(ssx);
    ssx_min=min(ssx);
    ssy=reshape(s(:,2),n,sM);
    ssy=max(ssy);
    ssy_max=max(ssy);
    ss=reshape(ones(n,1)*(ssx<=0&ssy>2*pi),sM*n,1);
    ss=logical(ss);
    
    ttx=reshape(t(:,1),n,tM);
    ttx=max(ttx);
    tty=reshape(t(:,2),n,tM);
    tty=min(tty);
    tt=reshape(ones(n,1)*(ttx>ssx_min+2*pi&tty<=ssy_max-2*pi),tM*n,1);
    tt=logical(tt);
    
    
    sN=round(sum(ss)/n); tN=round(sum(tt)/n);
    if sN>0&&tN>0
    edge=reshape((1:sN*n)',n,sN);
    edge=[(1:sN*n)',reshape(edge([2:end,1],:),sN*n,1)];
    in_t=inpoly(t(tt,:)+ones(tN*n,1)*[-2*pi 2*pi],s(ss,:),edge);
    in_t=reshape(in_t,n,tN);
    in_t=ones(n,1)*max(in_t);
    in_t=reshape(in_t,n*tN,1);
    in=zeros(n*tM,1);
    in(tt)=in_t;
    t=t(~in,:);
    end
    
        % Bottom Right Source, Top Left Target
    tM=round(length(t)/n);
    
    ssx=reshape(s(:,1),n,sM);
    ssx=max(ssx);
    ssx_max=max(ssx);
    ssy=reshape(s(:,2),n,sM);
    ssy=min(ssy);
    ssy_min=min(ssy);
    ss=reshape(ones(n,1)*(ssx>2*pi&ssy<=0),sM*n,1);
    ss=logical(ss);
    
    ttx=reshape(t(:,1),n,tM);
    ttx=min(ttx);
    tty=reshape(t(:,2),n,tM);
    tty=max(tty);
    tt=reshape(ones(n,1)*(ttx<=ssx_max-2*pi&tty>ssy_min+2*pi),tM*n,1);
    tt=logical(tt);
    
    
    sN=round(sum(ss)/n); tN=round(sum(tt)/n);
    if sN>0&&tN>0
    edge=reshape((1:sN*n)',n,sN);
    edge=[(1:sN*n)',reshape(edge([2:end,1],:),sN*n,1)];
    in_t=inpoly(t(tt,:)+ones(tN*n,1)*[2*pi -2*pi],s(ss,:),edge);
    in_t=reshape(in_t,n,tN);
    in_t=ones(n,1)*max(in_t);
    in_t=reshape(in_t,n*tN,1);
    in=zeros(n*tM,1);
    in(tt)=in_t;
    t=t(~in,:);
    end
    
end