function [ nV_tc] = OrthogonalizeTimeComponents(V_tc)

NbComp=size(V_tc,2);
nV_tc=V_tc;
for a=2:NbComp
    for b=1:(a-1)
         nV_tc(:,a)=nV_tc(:,a)-nV_tc(:,b)*nV_tc(:,b)'*V_tc(:,a)/norm(nV_tc(:,b))^2;
    end
end


end
