function nq = slerpInterpolation(oldSamp,oq,newSamp)

   oq = adjustQuatsForBlending(oq);

   nq = zeros(4,size(newSamp,2));

   for s = 1:numel(oldSamp)-1
      
      [~,c1]=find(newSamp>=oldSamp(s),1);
      [~,c2]=find(newSamp>=oldSamp(s+1),1);
      if isempty(c2)
         c2 = numel(newSamp);
      end
      
%       nqs = real(slerp(oq(:,s),oq(:,s+1),newSamp(c1:c2)-oldSamp(s),0.01));
      
      nqs = C_slerp(oq(:,s),oq(:,s+1),newSamp(c1:c2)-oldSamp(s));
      
      nq(:,c1:c2)=nqs;
      
   end

   nq = quatnormalize(nq);
   
end