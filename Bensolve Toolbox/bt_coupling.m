function [vals,B,Bs] = bt_coupling(Y,Ys,c,optdir)
	% -- function for internal use
	%
	%    see also http://tools.bensolve.org/files/manual.pdf

	%min: optdir=1
	%max: optdir=-1
	Y0=Y(1,:);
	Y1=Y(2:end-1,:);
	Yq=Y(end,:);
	Ys0=Ys(1,:);
	Ys1=Ys(2:end-1,:);
	Ysq=Ys(end,:);
	c1=c(1:end-1,1);
	cq=c(end,1);
	abscq=max(cq,-cq);
	vals=optdir*(cq*Y1'*Ys1 + Yq'*(cq/abscq*Ys0 - c1'*Ys1) - Y0'*Ysq);
	B=optdir*[Yq'*(cq/abscq),cq*Y1'-Yq'*c1',-Y0'];
	Bs=optdir*[-Ysq',cq*Ys1',cq/abscq*Ys0' - Ys1'*c1];
end

