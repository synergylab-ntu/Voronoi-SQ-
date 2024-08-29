function c_helptext(filename)
	n=strfind(filename,'.');
	fstr=filename(1:n-1);
	fstr=['./doc/',fstr,'_help.tex'];
	
	fid = fopen (filename);
	k=0;
	flag=0;
	out='';
	while(1)
		k=k+1;
		str = fgetl (fid);
		if str==-1
			break;
		end
		nn=strfind(str,'% **');
		if ~isempty(nn)
			out=[out,sprintf('\\subsubsection{%s}\n',str(nn(1)+5:end))]; %#ok
		end
		if contains(str,'% --')
			flag=1;
			out=[out,sprintf('\\begin{helptext}\n')];%#ok
			n=strfind(str,'%');
			if ~isempty(n)
				out=[out,sprintf('%s\n',str(n(1)+1:end))];%#ok
			end
		elseif flag==1
			n=strfind(str,'%');
			if ~isempty(n)
				if ~contains(str,'see also') && ~contains(str,'for further details') && ~contains(str,'For the')
					out=[out,sprintf('%s\n',str(n(1)+1:end))];%#ok
				else
					out=[out,sprintf('\\end{helptext}\n')];%#ok
					flag=0;
				end
			else
				out=[out,sprintf('\\end{helptext}\n')];%#ok
				flag=0;
			end
		end
	end
	fclose (fid);
	
	fid1=fopen(fstr,'w');
	fwrite(fid1,out);
	fclose (fid1);
end