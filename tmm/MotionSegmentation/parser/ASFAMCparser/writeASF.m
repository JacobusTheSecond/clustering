function writeASF(skel,filename)
% writeASF(skel,filename)

if (length(skel.nodes)<1)
    error('Can''t write empty skeleton.');
end

cl1 = clock;

%%%%%%%%%%%%%% applies to BVH-style input files only
if (strcmpi(skel.fileType,'BVH')) % convert BVH Euler angles into ASF-style Euler angles by reversing rotation order
    for k = skel.animated'
        skel.nodes(k).rotationOrder = fliplr(skel.nodes(k).rotationOrder);
    end
end

nodes = skel.nodes;
if ~strcmpi(nodes(1).boneName,'root')
    warning(['Bone 1 (' nodes(1).boneName ') renamed to "root"!']);
    nodes(1).boneName = 'root';
end

%%%%%%%%%%%%%% open ASF file
fid = fopen(filename,'wt');
if fid ~= -1
    %%%%%%% write comment
    fprintf(fid,'# ASF file generated by Tido''s ASF writer for MATLAB.\n');
    
    %%%%%%% write version
    if (isempty(skel.version))
        skel.version = '1.10';
    end
    fprintf(fid,':version %s\n',skel.version);
    
    %%%%%%% write name
    if (isempty(skel.name))
       skel.name = skel.fileType;
    end
    fprintf(fid,':name %s\n',skel.name);
    
    %%%%%%% write units
    fprintf(fid,':units\n  mass %.20g\n  length %.20g\n  angle %s\n',skel.massUnit, skel.lengthUnit, skel.angleUnit);
    
    %%%%%%% write documentation
    fprintf(fid,':documentation\n');
    for k=1:length(skel.documentation)
        fprintf(fid,'%s\n',skel.documentation{k,1});
    end

    %%%%%%% write root
    fprintf(fid,':root\n');
    fprintf(fid,'  axis %s\n',lower(nodes(1).rotationOrder));
    fprintf(fid,'  order ');
    for k=1:size(skel.nodes(1).DOF,1)
        fprintf(fid,'%s ',lower(nodes(1).DOF{k,1}));
    end
    fprintf(fid,'\n');
    fprintf(fid,'  position %.20g %.20g %.20g\n',nodes(1).offset  * skel.lengthUnit);
    fprintf(fid,'  orientation %.20g %.20g %.20g\n',skel.rootRotationalOffsetEuler);
        
    %%%%%%% write bonedata
    fprintf(fid,':bonedata\n');
    for k=2:length(nodes)
        fprintf(fid, '  begin\n');
        fprintf(fid, '    id %d\n',k-1);
        fprintf(fid, '    name %s\n',nodes(k).boneName);
        fprintf(fid, '    direction %.20g %.20g %.20g\n',nodes(k).direction);
        fprintf(fid, '    length %.20g\n',nodes(k).length * skel.lengthUnit);
        fprintf(fid, '    axis %.20g %.20g %.20g %s\n',lower(nodes(k).axis),lower(nodes(k).rotationOrder));
        if (size(nodes(k).DOF,1)>0)
            fprintf(fid,'    dof ');
            for j = 1:size(nodes(k).DOF,1)
                fprintf(fid, '%s ',lower(nodes(k).DOF{j,1}));
            end
            fprintf(fid,'\n');
        end
        if (size(nodes(k).limits,1)>0)
            fprintf(fid,'    limits ');
            for j = 1:size(nodes(k).limits,1)
                if j>1
                    fprintf(fid,'           ');
                end
                fprintf(fid, '(%.20g %.20g)\n',nodes(k).limits(j,:));
            end
        end        
        fprintf(fid, '  end\n');
    end
    
    %%%%%%% write hierarchy
    fprintf(fid,':hierarchy\n');
    fprintf(fid, 'begin\n');
    for k=1:length(nodes)
        if (length(nodes(k).children)>0)
            fprintf(fid,'  %s ',nodes(k).boneName);
            for j = 1:length(nodes(k).children)
                fprintf(fid,'%s ',nodes(nodes(k).children(j)).boneName);
            end
            fprintf(fid,'\n');
        end
    end
   fprintf(fid, 'end\n');
    
    %%%%%%% write skin
    if (size(skel.skin,1)>0)
        fprintf(fid,':skin ');
        for k=1:size(skel.skin,1)
            fprintf(fid,'%s\n',skel.skin{k});
        end
    end
    
    fclose(fid);        
else
    error('Could not open ASF file.');
end

t=etime(clock,cl1);
disp(['Wrote ASF file ' filename ' in ' num2str(t) ' seconds.']);


