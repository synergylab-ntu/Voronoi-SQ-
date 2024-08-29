function C = bt_polyh_delmult(A,tol)
	% -- function for internal use
	%
	%    see also http://tools.bensolve.org/files/manual.pdf
	
	[m,n]=size(A);
	B=zeros(m,n);
	k=0;
	for i=1:n
		if ~vectorismember(A(:,i),B(:,1:k),tol)
			k=k+1;
			B(:,k)=A(:,i);
		end
	end
	C=B(:,1:k);
end

function tf = vectorismember(x, S, tol)
	tf=0;
	for i=1:size(S,2)
		if norm(x-S(:,i)) < tol
			tf=1;
			break;
		end
	end
end