function data = readH5acc(filename,varargin)

dopts.DEBUG = true;

switch nargin
    case 1
        options = dopts;
    case 2 
        options = mergeOptions(varargin{1},dopts);
    otherwise
        error('Wrong num of args!');
end

info = h5info(filename);

nSensors = numel(info.Groups);

for i=1:nSensors
   fieldname = info.Groups(i).Name;
   fieldname = strrep(fieldname,'/','');
   fieldname = strrep(fieldname,'-','');
   
   data.sensor.(fieldname) = struct();
   
   nModalities = size(info.Groups(i).Groups(1).Datasets,1);
   
   for j=1:nModalities
        
        data.sensor.(fieldname).(info.Groups(i).Groups(1).Datasets(j).Name) ...
            = h5read(filename, ...
                [info.Groups(i).Groups(1).Name '/'  ...
                 info.Groups(i).Groups(1).Datasets(j).Name]);
   end
   
end

data.samplingRate = info.Groups(1).Attributes(1).Value;

fieldname = info.Groups(1).Name;
fieldname = strrep(fieldname,'/','');
fieldname = strrep(fieldname,'-','');

data.nframes = size( data.sensor.(fieldname).(info.Groups(1).Groups(1).Datasets(1).Name) , 2);

if options.DEBUG

    for j=1:nModalities
        
        cmod = info.Groups(i).Groups(1).Datasets(j).Name;
        
        figure();
        fields = fieldnames(data.sensor);
        for i=1:nSensors
            subplot(nSensors,1,i);
            plot(data.sensor.(fields{i}).(cmod)');
            grid on;
            title(['Sensor ' fields{i} ' ' cmod]);
        end
    end
    
end

end