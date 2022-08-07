invokeGapScript:=function(dir,name)
	local s;
	s:=ConcatenationString("labels_",name);
	s:=ConcatenationString(s,".g");
	s:=ConcatenationString(dir,s);
	Read(s);
end;

################################
####FINDING REGULAR ELEMENTS####
################################

isWellGenerated:=function(W)
	local r,n,h,i,deg,cod;
	r:=W.type[1];
	if r.series<>"ST" then
		return true;
	fi;
	deg:=ReflectionDegrees(W);
	n:=Length(deg);
	h:=Maximum(deg);
	Sort(deg);
	cod:=ReflectionCoDegrees(W);
	Sort(cod);
	for i in [1..n] do
		if deg[i]+cod[n+1-i]<>h then
			return false;
		fi;
	od;
	return true;
end;

allReflections:=function(W)
	#Christians Code
	local base,i,exp;

	base:=Set(Reflections(W));
	for i in [1..Length(base)] do 
		exp:=2;
		while not base[i]^exp=() do 
			Add(base,base[i]^exp); 
			exp:=exp+1;
		od;
	od;
	return base;
end;

candidate:=function(W,rMat)
	#Christians Code
	local var,h,i;

	Sort(ReflectionDegrees(W));
	h:=Reversed(ReflectionDegrees(W))[1];
	for i in [1..h] do
		var:=NullspaceMat(rMat-E(h)^i*IdentityMat(W.nbGeneratingReflections));
		if Gcd([i,h])=1 and not var=[] then 
			return(var[1]); 
		fi;
	od;
end;

isRegular:=function(W,rMat)
	#Christians Code
	local refls,cand,i;

	refls:=allReflections(W);
	cand:=candidate(W,rMat);
	for i in [1..Length(refls)] do 
		if cand-cand*MatXPerm(W,refls[i])=NullMat(1,W.nbGeneratingReflections)[1] then 
			return false;
		fi;
	od;
	return true;
end;

regularElement:=function(W)
	#Christians Code
	local rMat;

	#for avoiding conflicts with "Powermap not bound."
	ReflectionLength(W,());

	rMat:=Product(W.matgens);
	if isRegular(W,rMat) then
		return PermMatX(W,rMat);
	fi;
end;

###########################
####CREATING NC-LATTICE####
###########################

isSmaller:=function(W,u,v)
	if ReflectionLength(W,v)=(ReflectionLength(W,u)+ReflectionLength(W,u^-1*v)) then
		return true;
	else
		return false;
	fi;
end;

atoms:=function(W,reg)
	local refls,ats,r;

	refls:=allReflections(W);
	ats:=[];
	for r in refls do
		if isSmaller(W,r,reg)=true then
		Add(ats,r);
		fi;
	od;
	return ats;
end;

latticeLayers:=function(W,verbose)
	local reg,ats,h,i,layers,l,u,a,cand;

	reg:=regularElement(W);
	layers:=[];
	h:=ReflectionLength(W,reg);
	ats:=atoms(W,reg);
	l:=[reg];
	Add(layers,l);

	for i in [1..h-2] do
		l:=[];
		for u in layers[i] do
			for a in ats do
				cand:=a^-1*u;
				if isSmaller(W,cand,reg) then
					if verbose=true then
						PrintTo("*stdout*",".");
					fi;
					Add(l,cand);
				fi;
			od;
		od;
		Add(layers,Set(l));
	od;
	Add(layers,ats);
	Add(layers,[()]);
	#if verbose=true then
	#  PrintTo("*stdout*","\n");
	#fi;
	return layers;
end;

latticeElements:=function(W,verbose)
	local l,u,layers,elements;

	elements:=[];
	layers:=latticeLayers(W,verbose);
	for l in layers do
		for u in l do
			if verbose=true then
				PrintTo("*stdout*",".");
			fi;
			Add(elements,u);
		od;
	od;
	if verbose=true then
		PrintTo("*stdout*","\n");
	fi;
	return Set(elements);
end;

#############################
####OUTPUTTING NC-LATTICE####
#############################

intVal:=function(s)
	if s="1" then
		return 1;
	elif s="2" then
		return 2;
	elif s="3" then
		return 3;
	elif s="4" then
		return 4;
	elif s="5" then
		return 5;
	elif s="6" then
		return 6;
	elif s="7" then
		return 7;
	elif s="8" then
		return 8;
	elif s="9" then
		return 9;
	elif s="0" then
		return 0;
	else
		return -1;
	fi;
end;

toInt:=function(word)
	local i,n,sum;
	sum:=0;
	n:=Length(word);
	for i in [1..n] do
		if intVal(SubString(word,i,i))<>-1 then
			sum:=sum+10^(n-i)*intVal(SubString(word,i,i));
		else
			return -1;
		fi;
	od;
	return sum;
end;

rn:=function(l,el)
	local s,i,u,list;
  
	if IsGroup(l) then
		list:=latticeElements(l,false);
	else 
		list:=l;
	fi;

	i:=0;
	for u in list do
		i:=i+1;
		if u=el then
			s:=ConcatenationString("s",String(i));
			return s;
		fi;
	od;
	return ConcatenationString("s",String(-1));
end;

re:=function(l,s)
	local list;
  
	if IsGroup(l) then
		list:=latticeElements(l,false);
	else
		list:=l;
	fi;

	return list[toInt(SubString(s,2,Length(s)))];
end;

printLattice:=function(W,dest,verbose)
	local elements,u,v;
  
	if isWellGenerated(W) then
		elements:=latticeElements(W,verbose);
		PrintTo(dest,"NC(",W.name,")\n");
		for u in elements do
			for v in elements do
				if isSmaller(W,u,v) then
					AppendTo(dest,rn(elements,u),",",rn(elements,v),"\n");
				fi;
			od;
		od;
		AppendTo(dest,"##\n");
	else 
		PrintTo("*stdout*","This group is not well-generated.\n");
	fi;
end;

printGroupTable:=function(W,dest,verbose)
	local elements,u,v,uv;
  
	PrintTo(dest,"NC(",W.name,")\n");
	elements:=latticeElements(W,false);
	for u in elements do
		if verbose=true then
			PrintTo("*stdout*", ".");
		fi;
		for v in elements do
			uv:=rn(elements,u^-1*v);
			if uv<>"s-1" then
				AppendTo(dest,rn(elements,u),",",rn(elements,v),",",uv,"\n");
			fi;
		od;
	od;
	if verbose=true then
		PrintTo("*stdout*","\n");
	fi;
	AppendTo(dest,"##\n");
end;

#######################
####LABELING CHAINS####
#######################

maximalElement:=function(W,lattice)
	local max,u,kMax,kAct;

	kMax:=0;
	for u in lattice do
		kAct:=ReflectionLength(W,u);
		if kAct>kMax then
			kMax:=kAct;
			max:=u;
		fi;
	od;
	return max;
end;

minimalElement:=function(W,lattice)
	local min,u,kMin,kAct;

	kMin:=W.nbGeneratingReflections+1;
	for u in lattice do
		kAct:=ReflectionLength(W,u);
		if kAct<kMin then
			kMin:=kAct;
			min:=u;
		fi;
	od;
	return min;
end;

upperNeighbors:=function(W,nc,p)
	local u,cand,res,k;

	res:=[];
	k:=ReflectionLength(W,p);
	for u in nc do
		if ReflectionLength(W,u)=k+1 and isSmaller(W,p,u) then
			Add(res,u);
		fi;
	od;
	return res;
end;

maximalChains:=function(W,nc,max,min)
	local lu,i,u,c,chains,subchains;
	chains:=[];
	lu:=upperNeighbors(W,nc,min);
	if max in lu then
		return [[min,max]];
	fi;
	for i in [1..Length(lu)] do
		if isSmaller(W,lu[i],max) then
			subchains:=maximalChains(W,nc,max,lu[i]);
			for u in subchains do
				c:=[min];
				Append(c,u);
				Add(chains,c);
			od;
		fi;
	od;
	return chains;
end;

chains:=function(W,nc,interval)
	local min,max;
	min:=minimalElement(W,interval);
	max:=maximalElement(W,interval);
	return maximalChains(W,nc,max,min);
end;

printChain:=function(nc,c,dest)
	local u;
	for u in c do
		AppendTo(dest,rn(nc,u));
		if Position(c,u)<>Length(c) then
			AppendTo(dest,"----");
		fi;
	od;
	AppendTo(dest,"\n");
end;

labelChain:=function(nc,c,dest)
	local i;
	AppendTo(dest,"    ");
	for i in [1..Length(c)-1] do
		AppendTo(dest,rn(nc,c[i]^-1*c[i+1]));
		if i<>Length(c)-1 then
			AppendTo(dest,"    ");
		fi;
	od;
	AppendTo(dest,"\n");
end;

printLabeledChains:=function(W,interval,dest,verbose)
	local c,nc;

	nc:=latticeElements(W,false);
	for c in chains(W,nc,interval) do
		labelChain(nc,c,dest);
		printChain(nc,c,dest);
		AppendTo(dest,"\n");
		if verbose=true then
			PrintTo("*stdout*",".");
		fi;
	od;
	if verbose=true then
		PrintTo("*stdout*","\n");
	fi;
end;