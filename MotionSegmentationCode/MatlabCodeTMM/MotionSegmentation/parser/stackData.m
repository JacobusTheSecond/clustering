function dataStacked = stackData(data,offset,varargin)

if offset==0
    dataStacked = data;
else
    offset = sort(offset);
    nrRows = size(data,1);
    nrCols = size(data,2);
    dataStacked = zeros(nrRows*(numel(offset)),nrCols);

    if nargin == 2

        x = 0:nrCols-1;
        leftPufferSize  = abs(min(0,offset(1)));
        rightPufferSize =     max(0,offset(end));
        xx = -leftPufferSize:nrCols+rightPufferSize-1;

        data = interp1(x,data',xx,'linear','extrap')';
        if size(data,2)==1
            data=data';
        end

    %     data = double([nan(nrRows,-offset(1)),data,nan(nrRows,offset(end))]);

    %     data = double([prePadd, ...
    %                     data, ...
    %                    postPadd]);
        id = 1;
        for i=1:numel(offset)
            rws                 = ((i-1)*nrRows)+1:i*nrRows;
    %         dataStacked(rws,:)  = double(circshift(data,[0 -offset(i)]));
            dataStacked(rws,:)  = data(:,id:id+nrCols-1);
            if i<numel(offset)
                id = id + abs(offset(i+1)-offset(i));
            end
        end

    elseif nargin==3
        startIDs = [varargin{1}(:)' nrCols+1];
        for i=2:numel(startIDs)
            data_tmp = data(:,startIDs(i-1):startIDs(i)-1);
            data_tmpStacked = stackData(data_tmp,offset);
            dataStacked(:,startIDs(i-1):startIDs(i)-1) = data_tmpStacked;
        end
    else
        error('Wrong number of argins');
    end
end
