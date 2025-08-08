function path = relativePath(path)
  aux_pwd = pwd;
  aux_pos = strfind(aux_pwd,filesep);
  path_pwd = [];
  for k = 1:length(aux_pos)-1
    path_pwd{k} = aux_pwd((aux_pos(k)+1):(aux_pos(k+1)-1));
  end
  path_pwd{end+1} = aux_pwd((aux_pos(k+1)+1):end);

  aux_file = path;
  aux_pos = strfind(aux_file,filesep);  
  if aux_pos(end) == length(aux_file)
    aux_file(end) = [];
    aux_pos(end) = [];
  end  
  path_file = [];
  for k = 1:length(aux_pos)-1
    path_file{k} = aux_file((aux_pos(k)+1):(aux_pos(k+1)-1));
  end
  path_file{end+1} = aux_file((aux_pos(k+1)+1):end);

  on_loop = true;
  k = 0;
  while on_loop
    k = k + 1;
    if (k <= numel(path_pwd)) && (k <= numel(path_file))
      if ~strcmpi(path_pwd{k},path_file{k})
        on_loop = false;
      end  
    else
      on_loop = false;
    end
  end
  k = k-1;

  path = repmat(['..' '/'],1,length(path_pwd)-k);

  for k = (k+1):length(path_file)
    path = [path path_file{k} '/'];
  end
end

  