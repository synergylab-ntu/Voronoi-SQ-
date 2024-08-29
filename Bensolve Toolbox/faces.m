function  carr = faces(P,k)
	% -- carr = faces(P,k)   l-faces of P
	%
	%    Input:
	%      P: polyhedron (polyh object)
	%    Optional input:
	%      k: dimension of faces (number)
	%         default: return all proper faces
	%    Output:
	%      carr: cell array of faces
	%    
	%    see also:
	%
	%    for further details, see http://tools.bensolve.org/files/manual.pdf

	P=P.eval;
	if nargin==1
		gra=hasse(P);
		carr=cell(1,size(gra,1)-2);
		for i=1:size(gra,1)-2
			carr{i}=gra{i+1,3};
		end
	else
		gra=hasse(P,k,k>(P.dim-1)/2);
		carr=cell(1,size(gra,1));
		cnt=0;
		for i=1:size(gra,1)
			if gra{i,4}==k
				cnt=cnt+1;
				carr{cnt}=gra{i,3};
			end
		end
		carr=carr(1:cnt);
	end
end
