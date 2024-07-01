function [quats_out,q_ref,expmaps_out] = adjustQuatsForBlending(quats)

if iscell(quats)
    q_ref       = cell(numel(quats),1);
    quats_out   = quats;
    if nargout>=3
        expmaps_out = cell(numel(quats),1);
    end
    for j=1:numel(quats)

        if size(quats{j},1)~=4
            error('Input must be a cell array of 4xN quaternions or matrix of size (4M)xN!');
        end

        quats_tmp = quats{j};

        A = quats_tmp*quats_tmp';
        [V,D] = eig(A);

        q_ref{j}     = V(:,end);
        if(q_ref{j}(1)==-1)
           q_ref{j} = - q_ref{j};
        end
        
        q_ref_rep    = q_ref{j}(:,ones(size(quats_tmp,2),1));

        d   = sum(quats_tmp.*q_ref_rep);
        idx = acos(d)>pi/2;
        quats_tmp(:,idx) = -quats_tmp(:,idx);

        quats_out{j} = quats_tmp;

        if nargout>=3
            expmaps_out{j} = quatlog_refquat(quats_out{j},q_ref{j});
        end
    end
else
    if mod(size(quats,1),4)~=0
        error('Input must be a cell array of 4xN quaternions or matrix of size (4M)xN!');
    end
    q_ref = zeros(4,size(quats,1)/4);
    quats_out = zeros(size(quats));
    if nargout>=3
        expmaps_out = zeros(size(quats,1)/4*3,size(quats,2));
    end
    for j=1:size(quats,1)/4

        quats_tmp = quats(j*4-3:j*4,:);

        A = quats_tmp*quats_tmp';
        [V,D] = eig(A);

        q_ref(:,j) = V(:,end);
        
        if(q_ref(1,j)==-1)
           q_ref(:,j) = - q_ref(:,j);
        end
        
        q_ref_rep    = q_ref(:,j*ones(size(quats_tmp,2),1));

        d   = sum(quats_tmp.*q_ref_rep);
        idx = acos(d)>pi/2;
        quats_tmp(:,idx) = -quats_tmp(:,idx);

        quats_out(j*4-3:j*4,:) = quats_tmp;

        if nargout>=3
            expmaps_out(j*3-2:j*3,:) = quatlog_refquat(quats_out(j*4-3:j*4,:),q_ref(:,j));
        end
    end
end