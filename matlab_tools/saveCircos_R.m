function saveCircos_R(T,H,varNames,fname)

T = round(T*100)/100;
T(isnan(T)) = 0;

H = round(H*100)/100;
H(isnan(H)) = 0;

D = size(T,1);
assert(size(T,2)==D);
assert(length(H)==D);

% write transfer entropy

Tstring = cell(D+1,D);
for x = 1:D
 for y = 1:D
  Tstring{x+1,y} = num2str(abs(T(x,y)));
 end
end

Tstring(1,:) = varNames;

fid = fopen(strcat(fname(1:end-4),'_TE',fname(end-3:end)),'w');
format=repmat('%s ',[1,D]);
format = strcat(format,'\n');
for x = 1:D+1
 wstring = Tstring(x,:);
 fprintf(fid,format,wstring{1:end});
end
fclose(fid);

% write entropy

hname = strcat(fname(1:end-4),'_H',fname(end-3:end));
csvwrite(hname,H);

