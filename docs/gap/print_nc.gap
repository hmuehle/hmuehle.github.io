##############################
####NONCROSSING PARTITIONS####
##############################

is_well_generated:=function(W)
	# input: W .. complex reflection group
	# output: true if W is well-generated
	local r,n,h,i,deg,cod;
	r := W.type[1];
	if r.series<>"ST" then
		return true;
	fi;
	n := Length(r.degrees);
	h := Maximum(r.degrees);
	deg := ReflectionDegrees(W);
	Sort(deg);
	cod := ReflectionCoDegrees(W);
	Sort(cod);
	for i in [1..n] do
		if deg[i] + cod[n+1-i] <> h then
			return false;
		fi;
	od;
	return true;
end;

all_reflections:=function(W)
	# input: W .. complex reflection group
	# output: set of all reflections of W
	local base,i,exp;

	base := Set(Reflections(W));
	for i in [1..Length(base)] do 
		exp := 2;
		while not base[i]^exp = () do 
			Add(base,base[i]^exp); 
			exp := exp+1;
		od;
	od;
	return base;
end;

candidate:=function(W,rMat)
	# input: W .. complex reflection group, rMat .. matrix of some group element
	# output: an eigenvector of rMat
	local var,h,i;

	Sort(ReflectionDegrees(W));
	h := Reversed(ReflectionDegrees(W))[1];
	for i in [1..h] do
		var := NullspaceMat(rMat-E(h)^i*IdentityMat(W.nbGeneratingReflections));
		if Gcd([i,h]) = 1 and not var = [] then 
			return(var[1]); 
		fi;
	od;
end;

is_regular:=function(W,rMat)
	# input: W .. complex reflection group, rMat .. matrix of some group element
	# output: true, if rMat corresponds to a regular element
	local refls,cand,i;

	refls := all_reflections(W);
	cand := candidate(W,rMat);
	for i in [1..Length(refls)] do 
		if cand - cand * MatXPerm(W,refls[i]) = NullMat(1,W.nbGeneratingReflections)[1] then 
			return false;
		fi;
	od;
	return true;
end;

regular_element:=function(W)
	# input: W .. complex reflection group
	# output: a regular element of W, if possible
	local rMat;

	#for avoiding conflicts with "Powermap not bound."
	ReflectionLength(W,());

	rMat := Product(W.matgens);
	if is_regular(W,rMat) then
		return PermMatX(W,rMat);
	fi;
end;

is_abs_leq:=function(W,u,v)
	# input: W .. complex reflection group, u,v .. elements of W
	# output: true, if u is smaller or equal to v in the absolute order of W
	return ReflectionLength(W,v) = (ReflectionLength(W,u) + ReflectionLength(W,u^-1*v));
end;

atoms:=function(W,reg)
	# input: W .. complex reflection group, reg .. regular element of W
	# output: reflections of W below reg in absolute order
	local refls,ats,r;

	refls := all_reflections(W);
	ats := [];
	for r in refls do
		if is_abs_leq(W,r,reg) = true then
		Add(ats,r);
		fi;
	od;
	return ats;
end;

lattice_layers:=function(W,verbose)
	# input: W .. complex reflection group, verbose .. switch to toggle output
	# output: ranks of the noncrossing partition lattice of W
	local reg,ats,h,i,layers,l,u,a,cand,run,s;
	
	if verbose then
		Exec("date");
		PrintTo("*stdout*","Computing Lattice Layers.");	
	fi;
	reg := regular_element(W);
	layers := [];
	h := ReflectionLength(W,reg);
	ats := atoms(W,reg);
	l := [reg];
	Add(layers,l);
	run := 0;

	for i in [1..h-2] do
		if verbose then
			PrintTo("*stdout*","\n  Computing Layer ",String(h-i),"\n");
			PrintTo("*stdout*","\tTotal: ",String(Length(layers[i])*Length(ats)),"\n");
		fi;
		l := [];
		for u in layers[i] do
			for a in ats do
				cand := a^-1*u;
				if is_abs_leq(W,cand,reg) then
					Add(l,cand);
				fi;
				if verbose then
					run := run + 1;
					PrintTo("*stdout*","\r\tDone:  ",String(run));
				fi;
			od;
		od;
		Add(layers,Set(l));
	od;
	Add(layers,ats);
	Add(layers,[()]);
	if verbose then
		PrintTo("*stdout*","\n");
	fi;
	return layers;
end;

lattice_elements:=function(W,verbose)
	# input: W .. complex reflection group, verbose .. switch to toggle output
	# output: elements of the noncrossing partition lattice of W
	local l,u,layers,elements,run;

	elements := [];
	layers := lattice_layers(W,verbose);
	run := 0;
	for l in layers do
		for u in l do
			run := run + 1;
			Add(elements,u);
		od;
	od;
	return Set(elements);
end;

############################################
####ELEMENT CONVERSION (PRETTY PRINTING)####
############################################

int_val:=function(s)
	# input: s .. string
	# output converts s to an integer, if s is a digit
	if s = "1" then
		return 1;
	elif s = "2" then
		return 2;
	elif s = "3" then
		return 3;
	elif s = "4" then
		return 4;
	elif s = "5" then
		return 5;
	elif s = "6" then
		return 6;
	elif s = "7" then
		return 7;
	elif s = "8" then
		return 8;
	elif s = "9" then
		return 9;
	elif s = "0" then
		return 0;
	else
		return -1;
	fi;
end;

to_int:=function(word)
	# input: word .. string
	# output: converts word to an integer, if possible
	local i,n,sum;
	sum := 0;
	n := Length(word);
	for i in [1..n] do
		if int_val(SubString(word,i,i)) <> -1 then
			sum := sum + 10^(n-i) * int_val(SubString(word,i,i));
		else
			return -1;
		fi;
	od;
	return sum;
end;

rn:=function(l,el)
	# input: l .. list of elements, el .. element of l
	# output: "standard" name of el with respect to l
	local s,i,u,list;
  
	if IsGroup(l) then
		list := lattice_elements(l,false);
	else 
		list := l;
	fi;

	i := 0;
	for u in list do
		i := i + 1;
		if u = el then
			s := ConcatenationString("s",String(i));
			return s;
		fi;
	od;
	return ConcatenationString("s",String(-1));
end;

re:=function(l,s)
	# input: l .. list of elements, s .. string
	# output: element in l corresponding to s
	local list;
  
	if IsGroup(l) then
		list := lattice_elements(l,false);
	else
		list := l;
	fi;

	return list[to_int(SubString(s,2,Length(s)))];
end;

print_lattice_from_index:=function(W,ind,dest,verbose)
	# input: W .. complex reflection group, ind .. index of element to start with dest .. path, verbose .. switch to toggle output
	# output: produces a file at destination dest that contains the cover relations of the noncrossing partition lattice of W
	local elements,i,u,v,run;

	if is_well_generated(W) then
		elements := lattice_elements(W,verbose);
		if verbose then
			Exec("date");
			PrintTo("*stdout*","Printing Lattice.\n");
		fi;
		run := 0;
		if verbose then
			PrintTo("*stdout*","\tTotal: ",String(Length(elements)*(Length(elements)-ind+1)),"\n");
		fi;
		for i in [ind..Length(elements)] do
			u := elements[i];
			for v in elements do
				if is_abs_leq(W,u,v) and ReflectionLength(W,v) = ReflectionLength(W,u) + 1 then
					AppendTo(dest,rn(elements,u),",",rn(elements,v),"\n");
				fi;
				if verbose then
					run := run + 1;
					PrintTo("*stdout*","\r\tDone:  ",String(run));
				fi;
			od;
		od;
		AppendTo(dest,"##\n");
	else 
		PrintTo("*stdout*","This group is not well-generated.\n");
	fi;
	if verbose then
		PrintTo("*stdout*","\n");
	fi;
end;

print_lattice:=function(W,dest,verbose)
	local elements,i,j,u,v;
	if is_well_generated(W) then
		PrintTo(dest,"NC(",W.name,")\n");
		print_lattice_from_index(W,1,dest,verbose);
	else 
		PrintTo("*stdout*","This group is not well-generated.\n");
	fi;
end;
