function [VSol,QSol,sSol,dSol,finalAlpha,alphas,diffRec,t,stateCurve]=hemMachinePFmultiStage(SimData,SysData,SysPara,x,pShare)
% Multi-stage computation of power flow (extended including synchronous generators)
%
% FUNCTION hemMachinePFmultiStage
%
% Author: Rui Yao <ruiyao@ieee.org>
%
% Copyright (C) 2021, UChicago Argonne, LLC. All rights reserved.
% 
% 
global TaskID TaskNum

setting_func_hemMachinePFmultiStage;

[bus,sw,pv,pq,shunt,line,ind,zip,syn,exc,tg,agc,cac,cluster]=unfoldSysData(SysData);
[nState,idxs]...
    =getIndexDyn(SysData);
[V0,Q0,s0,d0,w,eq1,eq2,ed1,ed2,psid,psiq,Pm0,Ef0,Vavrm,Vavrr,Vavrf,Vavrref,tgovg,tgovm,tgovmech,f,dpg,qplt,vg]=unfoldX(x,SysData);
[pqIncr,pvIncr,Rind0,Rind1,Reind0,Reind1,Rzip0,Rzip1,Ytr,~,Ysh0,Ysh1,VspSq2,~,~,~,~,Tmech1,Varref1,Ef1,Pm1,Eq11]=unfoldSysPara(SysPara);
[maxAlpha,segAlpha,dAlpha,nlvl,taylorN,alphaTol,diffTol,diffTolMax,method]=unfoldSimData(SimData);

if isempty(dAlpha)
    dAlpha=0.1;
    if dAlpha>maxAlpha/10
        dAlpha=maxAlpha/10;
    end
end

if isfield(SysPara,'nIslands')&&isfield(SysPara,'islands')&&isfield(SysPara,'refs')
    nIslands=SysPara.nIslands;islands=SysPara.islands;refs=SysPara.refs;
else
    [nIslands,islands,refs]=searchIslands(bus(:,1),line(:,1:2));
end

nbus=size(bus,1);
nline=size(line,1);

pqx=pq;
pvx=pv;
shuntx=shunt;
indx=ind;
zipx=zip;

V0x=V0;
s0x=s0;
Q0x=Q0;
d0x=d0;
Ef0x=Ef0;
Rind0x=Rind0;
Reind0x=Reind0;
Rzip0x=Rzip0;
Ysh0x=Ysh0;
Pm0x=Pm0;
Vsp0x=real(V0.*conj(V0));

alphaConfirm=0;

busType=zeros(nbus,1);
if ~isempty(pv);busType(pv(:,1))=1;end
if ~isempty(sw);busType(sw(:,1))=2;end
Vsp=zeros(nbus,1);
Vsp(busType==2)=sw(:,4).*exp(1j*sw(:,5))-V0(busType==2,1);

vc=zeros(nbus,1);
dc=zeros(size(syn,1),1);
qc=zeros(size(Q0,1),1);
Jc=zeros(size(syn,1),1);
Kc=zeros(size(syn,1),1);
sc=zeros(size(ind,1),1);
dsc=zeros(size(ind,1),1);
dyn=zeros(size(ind,1),1);

alphas=0;
noMove=0;

if exist('TaskID','var')&&exist('TaskNum','var')&&~isempty(TaskID)&&~isempty(TaskNum)
taskTag=num2str(TaskID(1));
else
taskTag='';
end

iden=[datestr(clock,30),'_',taskTag,'K',num2str(randi([100000,899999])),'_init'];
jac=getJacobianMod(nbus,nline,bus,sw,pvx,pqx,shuntx,line,zipx,indx,s0,V0,ones(size(ind,1),1));
initJacCond=condest(jac);

t=zeros(1,0);
stateCurve=zeros(nState,0);

while alphaConfirm<1-alphaTol/1000
    alpha=min([1-alphaConfirm,segAlpha]);
    alphax=alpha/alphaPreMult;
    
    x=zeros(nState,1);
    x(idxs.vIdx)=V0x;x(idxs.sIdx)=s0x;x(idxs.deltaIdx)=d0x;x(idxs.qIdx)=Q0x;x(idxs.efIdx)=Ef0x;x(idxs.pgIdx)=Pm0x;
    SysDatax=foldSysData(bus,sw,pvx,pqx,shuntx,line,indx,zipx,syn,exc,tg,agc,cac,cluster);
    SysParax=foldSysPara(pqIncr,pvIncr,Rind0x,Rind1,Reind0x,Reind1,Rzip0x,Rzip1,Ytr,[],Ysh0x,Ysh1,[Vsp0x,VspSq2(:,2)],[],[],[],[],[],[],Ef1,Pm1);
    SysParax.iden=iden;
    if exist([iden,'.mat'],'file')
        if ~exist('p_amd','var')
            load([iden,'.mat']);
        end
        SysParax.p_amd=p_amd;
    else
        SysParax.p_amd=[];
    end
    SysParax.nIslands=nIslands;
    SysParax.islands=islands;
    SysParax.refs=refs;   
    
    [V,W,Q,s,d,JG,KG]=hemMachinePFSalientcontinueSep(SimData,SysDatax,SysParax,x,pShare,Vsp);
    
    [Va3,Vb3]=getPadeApproxHelm(V,0);
    [da3,db3]=getPadeApproxHelm(d,0);
    [Qa3,Qb3]=getPadeApproxHelm(Q,find(busType==0));
    [Sa3,Sb3]=getPadeApproxHelm(s,0);
    ds=s(:,2:end).*repmat(1:(size(s,2)-1),size(s,1),1);
    [ds3a,ds3b]=getPadeApproxHelm(ds,0);
    idxNanV=find(isnan(Va3(:,1))|isnan(Vb3(:,1)));Va3(idxNanV,:)=0;Vb3(idxNanV,:)=0;Va3(idxNanV,1:3)=V(idxNanV,1:3);
    idxNanQ=find(isnan(Qa3(:,1))|isnan(Qb3(:,1)));Qa3(idxNanQ,:)=0;Qb3(idxNanQ,:)=0;Qa3(idxNanQ,1:3)=Q(idxNanQ,1:3);
    
    if ~isempty(ind)
        idxNanS=find(isnan(Sa3(:,1))|isnan(Sb3(:,1)));Sa3(idxNanS,:)=0;Sb3(idxNanS,:)=0;Sa3(idxNanS,1:2)=s(idxNanS,1:2);
        idxNandS=find(isnan(ds3a(:,1))|isnan(ds3b(:,1)));ds3a(idxNandS,:)=0;ds3b(idxNandS,:)=0;ds3a(idxNandS,1:2)=ds(idxNandS,1:2);
    end    
    
    if ~isempty(syn)
        idxNand=find(isnan(da3(:,1))|isnan(db3(:,1)));da3(idxNand,:)=0;db3(idxNand,:)=0;da3(idxNand,1:2)=d(idxNand,1:2);
    end
    %Temp: benchmarking with Cxx
%     CC=real(V);
%     DD=imag(V);
%     [CCa3,CCb3]=getPadeApproxHelm(CC,0);
%     [DDa3,DDb3]=getPadeApproxHelm(DD,0);
%     
%     X=[CC;DD;Q;s;d];
%     [Xa3,Xb3]=getPadeApproxHelm(X,0);
%     Xa=[CCa3;DDa3;Qa3;Sa3;da3];
%     Xb=[CCb3;DDb3;Qb3;Sb3;db3];
    
    %
    
    alphaLeft=0;
    alphaRight=alpha/alphaPreMult;
%     alphax=0.2;
    while 1
        
        vc=polyvalVec(Va3(:,end:-1:1),alphax)./polyvalVec([Vb3(:,end:-1:1),ones(size(Vb3,1),1)],alphax);
        
%         CCc=polyvalVec(CCa3(:,end:-1:1),alphax)./polyvalVec([CCb3(:,end:-1:1),ones(size(CCb3,1),1)],alphax);
%         CCc(isnan(CCc))=CC(isnan(CCc),1)+alphax*CC(isnan(CCc),2);
%         
%         DDc=polyvalVec(DDa3(:,end:-1:1),alphax)./polyvalVec([DDb3(:,end:-1:1),ones(size(DDb3,1),1)],alphax);
%         DDc(isnan(DDc))=DD(isnan(DDc),1)+alphax*DD(isnan(DDc),2);
%         vc=CCc+1j*DDc;
        
        dc=polyvalVec(da3(:,end:-1:1),alphax)./polyvalVec([db3(:,end:-1:1),ones(size(db3,1),1)],alphax);
        
        qc=polyvalVec(Qa3(:,end:-1:1),alphax)./polyvalVec([Qb3(:,end:-1:1),ones(size(Qb3,1),1)],alphax);
        qc(busType==0,:)=0;
                
        sc=polyvalVec(Sa3(:,end:-1:1),alphax)./polyvalVec([Sb3(:,end:-1:1),ones(size(Sb3,1),1)],alphax);
        dsc=polyvalVec(ds3a(:,end:-1:1),alphax)./polyvalVec([ds3b(:,end:-1:1),ones(size(ds3b,1),1)],alphax);
                
        Y=Ytr+sparse(1:nbus,1:nbus,Ysh0x+Ysh1*alphax,nbus,nbus);
        Ef=Ef0x+alphax*Ef1;
        Vsp2=Vsp0x+alphax*VspSq2(:,2);
        
        pqxx=pqx;
        pqxx(:,[4,5])=pqx(:,[4,5])+alphax*pqIncr;
        pvxx=pvx;
        if ~isempty(pvx);pvxx(:,4)=pvx(:,4)+alphax*pvIncr;end
        
        zipxx=zipx;
        if ~isempty(zipxx)            
            zipxx(:,5:10)=repmat(Rzip0x+alphax*Rzip1,1,6).*zipx(:,5:10);
        end
        
        indxx=indx;
        if ~isempty(indx)
            indxx(:,13)=indx(:,13)./(Reind0x+alphax*Reind1);
            indxx(:,15)=indx(:,15).*(Rind0x+alphax*Rind1);
            indxx(:,16)=indx(:,16).*(Rind0x+alphax*Rind1);
            indxx(:,17)=indx(:,17).*(Rind0x+alphax*Rind1);
            
            dsc=dsc.*dyn;
        end
        
        x=zeros(nState,1);
        x(idxs.vIdx)=vc;x(idxs.sIdx)=sc;x(idxs.deltaIdx)=dc;x(idxs.qIdx)=qc;x(idxs.efIdx)=Ef;
        dx=zeros(nState,1);
        dx(idxs.sIdx)=dsc;
        SysDatax=foldSysData(bus,sw,pvxx,pqxx,shunt,line,indxx,zipxx,syn,exc,tg,agc,cac,cluster);
        SysParax=foldSysPara([],[],[],[],[],[],[],[],Ytr,[],Ysh0x+Ysh1*alphax,[],[Vsp2,zeros(size(Vsp2))],[],[],[],[]);
        diff=checkEquationBalanceSyn(SysDatax,SysParax,x,dx,dyn);
        
        if max(abs(diff))<diffTol&&abs(alphax-alpha/alphaPreMult)<alphaTol/1000
            
            alphaLeft=alphaRight;
        else
            if max(abs(diff))<diffTol
                           
                alphaLeft=alphax;
                alphax=(alphaLeft+alphaRight)/2;
            else
                alphaRight=alphax;
                alphax=(alphaLeft+alphaRight)/2;
            end
        end
        
        if alphaRight-alphaLeft<alphaTol
            break;
        end
    end    
    
    diffParadigm=max(abs(diff))+1e-12;
    countCheckAlpha=1;
    alphax=alphaLeft*alphaPreMult;
    if alpha-alphax<alphaTol/10
        alphax=alpha;
    end
%     maxCount=10;
    while countCheckAlpha<maxCount
                
        vc=polyvalVec(Va3(:,end:-1:1),alphax)./polyvalVec([Vb3(:,end:-1:1),ones(size(Vb3,1),1)],alphax);
        
%         CCc=polyvalVec(CCa3(:,end:-1:1),alphax)./polyvalVec([CCb3(:,end:-1:1),ones(size(CCb3,1),1)],alphax);
%         CCc(isnan(CCc))=CC(isnan(CCc),1)+alphax*CC(isnan(CCc),2);
%         
%         DDc=polyvalVec(DDa3(:,end:-1:1),alphax)./polyvalVec([DDb3(:,end:-1:1),ones(size(DDb3,1),1)],alphax);
%         DDc(isnan(DDc))=DD(isnan(DDc),1)+alphax*DD(isnan(DDc),2);
%         vc=CCc+1j*DDc;
        
        dc=polyvalVec(da3(:,end:-1:1),alphax)./polyvalVec([db3(:,end:-1:1),ones(size(db3,1),1)],alphax);
        
        qc=polyvalVec(Qa3(:,end:-1:1),alphax)./polyvalVec([Qb3(:,end:-1:1),ones(size(Qb3,1),1)],alphax);
        qc(busType==0,:)=0;
        
        sc=polyvalVec(Sa3(:,end:-1:1),alphax)./polyvalVec([Sb3(:,end:-1:1),ones(size(Sb3,1),1)],alphax);
        dsc=polyvalVec(ds3a(:,end:-1:1),alphax)./polyvalVec([ds3b(:,end:-1:1),ones(size(ds3b,1),1)],alphax);
                
        Ef=Ef0x+alphax*Ef1;
        Vsp2=Vsp0x+alphax*VspSq2(:,2);
        
        pqxx=pqx;
        pqxx(:,[4,5])=pqx(:,[4,5])+alphax*pqIncr;
        pvxx=pvx;
        if ~isempty(pvx);pvxx(:,4)=pvx(:,4)+alphax*pvIncr;end
        
        zipxx=zipx;
        if ~isempty(zipxx)            
            zipxx(:,5:10)=repmat(Rzip0x+alphax*Rzip1,1,6).*zipx(:,5:10);
        end
        
        indxx=indx;
        if ~isempty(indx)
            indxx(:,13)=indx(:,13)./(Reind0x+alphax*Reind1);
            indxx(:,15)=indx(:,15).*(Rind0x+alphax*Rind1);
            indxx(:,16)=indx(:,16).*(Rind0x+alphax*Rind1);
            indxx(:,17)=indx(:,17).*(Rind0x+alphax*Rind1);
            
            dsc=dsc.*dyn;
        end
        
        x=zeros(nState,1);
        x(idxs.vIdx)=vc;x(idxs.sIdx)=sc;x(idxs.deltaIdx)=dc;x(idxs.qIdx)=qc;x(idxs.efIdx)=Ef;
        dx=zeros(nState,1);
        dx(idxs.sIdx)=dsc;
        SysDatax=foldSysData(bus,sw,pvxx,pqxx,shunt,line,indxx,zipxx,syn,exc,tg,agc,cac,cluster);
        SysParax=foldSysPara([],[],[],[],[],[],[],[],Ytr,[],Ysh0x+Ysh1*alphax,[],[Vsp2,zeros(size(Vsp2))],[],[],[],[]);
        diff=checkEquationBalanceSyn(SysDatax,SysParax,x,dx,dyn);
                
        if max(abs(diff))<1.05*diffParadigm
            break;
        end
        alphax=alphax*alphaPreMult;
        countCheckAlpha=countCheckAlpha+1;
    end
    if countCheckAlpha>=maxCount;alphax=alphax/alphaPreMult;end
    alpha=alphax;
    pqx(:,[4,5])=pqx(:,[4,5])+alpha*pqIncr;
    if ~isempty(pvx);pvx(:,4)=pvx(:,4)+alpha*pvIncr;end
    Rind0x=Rind0x+alpha*Rind1;
    Reind0x=Reind0x+alpha*Reind1;
    Rzip0x=Rzip0x+alpha*Rzip1;
    Ef0x=Ef0x+alpha*Ef1;
    Ysh0x=Ysh0x+alpha*Ysh1;
    Vsp0x=Vsp0x+alpha*VspSq2(:,2);
    
    if ~isempty(syn)
        nSyn=size(syn,1);
        synIdx=syn(:,1);
        Efd=zeros(nSyn,1);
        Efq=Ef0x;
        Rs=syn(:,7);
        Xd=syn(:,8);
        Xq=syn(:,13);
        cosd=cos(dc);
        sind=sin(dc);
        Cg=real(vc(synIdx));
        Dg=imag(vc(synIdx));
        Vd=sind.*Cg-cosd.*Dg;
        Vq=cosd.*Cg+sind.*Dg;
        Id=(Rs.*(Efd-Vd)+Xq.*(Efq-Vq))./(Rs.*Rs+Xq.*Xd);
        Iq=(-Xd.*(Efd-Vd)+Rs.*(Efq-Vq))./(Rs.*Rs+Xq.*Xd);
        Pm0x=Efq.*Iq+Efd.*Id;
    end
    
    alphaConfirm=alphaConfirm+alpha;
    addLog(['Step=',num2str(alphaConfirm),', added=',num2str(alpha),', (maxdiff<',num2str(diffTol),').'],'INFO');
    
    if alpha==0
        addLog('Step did not move!','INFO');
        noMove=noMove+1;
        if noMove>=MaxNoMove 
            zipxx=zipx;
            if ~isempty(zipxx)
                zipxx(:,5:10)=repmat(Rzip0x,1,6).*zipx(:,5:10);
            end
            [~,Ytrx,Yshx]=getYMatrix(nbus,line);
            [jac,jPIdx,jQIdx,jVIdx,jSIdx,busType]=getJacobianMod(nbus,nline,bus,sw,pvx,pqx,shuntx,line,zipxx,indx,sc,vc,zeros(size(ind,1),1),Ytrx,Yshx*alphaConfirm);
            jacCond=condest(jac);
            
            if jacCond>50*initJacCond
                [Vx,Dx]=eigs(jac,20,'sm');
                dDx=diag(Dx);
                [~,iSmall]=min(abs(dDx));
                part=Vx(:,iSmall).*conj(Vx(:,iSmall));
                partBus=zeros(size(bus,1),1);
                partBus(busType~=2)=partBus(busType~=2)+part(jPIdx);
                partBus(busType==0)=partBus(busType==0)+part(jQIdx);
                addLog('Singular likely!','INFO');
            end
            addLog('Reached consecutive max no move, exit!','INFO');
            break;
        end
        if diffTol>=diffTolMax
            addLog('Max DiffTol reached and not move, exit!','INFO');
            break;
        end
        diffTol=max([diffTol*diffMul,max(abs(diff))]);
        if diffTol>diffTolMax
            diffTol=diffTolMax;
        end
        addLog(['Enlarge tol! (Tol=',num2str(diffTol),')'],'INFO');
    else
        noMove=0;
        
        V0x=vc;
        d0x=dc;
        Q0x=qc;
        s0x=sc;
                    
        alphas=[alphas;alpha];
        
        dAlphax=dAlpha;
%         if dAlphax>alpha/10
%             dAlphax=alpha/10;
%         end
        
        alphax=0:dAlphax:alpha;
        if alpha-alphax(end)>alphaTol
            alphax=[alphax,alpha];
        end
        nt=size(alphax,2);
        
        vcx=polyvalVec(Va3(:,end:-1:1),alphax)./polyvalVec([Vb3(:,end:-1:1),ones(size(Vb3,1),1)],alphax);
        
        dcx=polyvalVec(da3(:,end:-1:1),alphax)./polyvalVec([db3(:,end:-1:1),ones(size(db3,1),1)],alphax);
        
        qcx=polyvalVec(Qa3(:,end:-1:1),alphax)./polyvalVec([Qb3(:,end:-1:1),ones(size(Qb3,1),1)],alphax);
        qcx(busType==0,:)=0;
        
        scx=polyvalVec(Sa3(:,end:-1:1),alphax)./polyvalVec([Sb3(:,end:-1:1),ones(size(Sb3,1),1)],alphax);
        
        if ~isempty(syn);efcx=repmat(Ef0x-alpha*Ef1,1,nt)+Ef1*alphax;else efcx=zeros(0,nt);end        
        
        xx=zeros(nState,nt);
        xx(idxs.vIdx,:)=vcx;xx(idxs.sIdx,:)=scx;xx(idxs.deltaIdx,:)=dcx;xx(idxs.qIdx,:)=qcx;xx(idxs.efIdx,:)=efcx;
        t=[t,alphax+alphaConfirm-alpha];
        stateCurve=[stateCurve,xx];
    end
end
addLog(['Total Step=',num2str(alphaConfirm),'.'],'INFO');

diffRec=max(abs(diff));
VSol=V0x;
QSol=Q0x;
sSol=s0x;
dSol=d0x;
finalAlpha=alphaConfirm;

if exist([iden,'.mat'],'file')
    delete([iden,'.mat']);
end
end
