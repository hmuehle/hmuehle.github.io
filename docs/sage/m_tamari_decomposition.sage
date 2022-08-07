#######################
#####M-TAMARI##########
#######################

def mTamari(m,n):
	"""
	Returns the m-Tamari lattice of parameter n.
	
	EXAMPLES::
	
		sage: mTamari(2, 3)
		Finite lattice containing 12 elements
	"""
	import sage.combinat.tamari_lattices as tam
	return tam.GeneralizedTamariLattice(m*n+1,n,m)
	
#####################
#####DECOMPOSING#####
#####################

def convert_to_height(p):
	"""
	Converts an m-Dyck path to its height sequence.
	
	The *height sequence* of an m-Dyck path `p` is the sequence of length `mn`, where the `i`-th entry is the y-coordinate of `p` at `x=i`.  This function converts an m-Dyck path to its height sequence.
	
	:param p: an m-Dyck path.
	
	EXAMPLES::
	
		sage: [convert_to_height(x) for x in mTamari(2,3)]
		[(1, 1, 2, 2, 3, 3),
		 (1, 1, 2, 3, 3, 3),
		 (1, 1, 3, 3, 3, 3),
		 (1, 2, 2, 2, 3, 3),
		 (1, 2, 2, 3, 3, 3),
		 (1, 2, 3, 3, 3, 3),
		 (1, 3, 3, 3, 3, 3),
		 (2, 2, 2, 2, 3, 3),
		 (2, 2, 2, 3, 3, 3),
		 (2, 2, 3, 3, 3, 3),
		 (2, 3, 3, 3, 3, 3),
		 (3, 3, 3, 3, 3, 3)]
	"""
	H = tuple([])
	h = 0
	for i in range(len(p)-1):
		if p[i]==1:
			h = h + 1
		if p[i]==0:
			H = H + tuple([h])
	return H

def convert_to_height_all(m,n):
	"""
	Converts all m-Dyck paths of parameter n to their corresponding height sequences.
	
	:param n: an integer.
	:param m: an integer.
	
	EXAMPLES::
	
		sage: convert_to_height_all(2,3)
		[(1, 1, 2, 2, 3, 3),
		 (1, 1, 2, 3, 3, 3),
		 (1, 1, 3, 3, 3, 3),
		 (1, 2, 2, 2, 3, 3),
		 (1, 2, 2, 3, 3, 3),
		 (1, 2, 3, 3, 3, 3),
		 (1, 3, 3, 3, 3, 3),
		 (2, 2, 2, 2, 3, 3),
		 (2, 2, 2, 3, 3, 3),
		 (2, 2, 3, 3, 3, 3),
		 (2, 3, 3, 3, 3, 3),
		 (3, 3, 3, 3, 3, 3)]
	"""
	return [convert_to_height(p) for p in mTamari(m,n)]

def convert_from_height(h,m,n):
	"""
	Converts a height sequence to the corresponding m-Dyck path of parameter n.
	
	:params h: a height sequence.
	:params m: an integer.
	:params n: an integer.
	
	EXAMPLES::
	
		sage: [convert_from_height(convert_to_height(p),2,3)==p for p in mTamari(2,3)]
		[True, True, True, True, True, True, True, True, True, True, True, True]
	"""
	p = tuple([])
	#initial up-steps
	for k in range(h[0]):
		p = p + tuple([1])
	for i in range(len(h)-1):
		p = p + tuple([0])
		for k in range(h[i+1]-h[i]):
			p = p + tuple([1])
	#final right-steps
	f = (m+1)*n-len(p)
	for i in range(f+1):
		p = p + tuple([0])
	return p

def decompose(h,m,n):
	"""
	Produces the strip-decomposition of an m-Dyck path as defined in Definition~3.9 of [KM15].
	
	Let `h` be the height sequence of an m-Dyck path `p`.  The *strip-decomposition* of `p` is the m-tuple (p_1,p_2,\ldots,p_m) of Dyck paths, where the height sequence of `p_i` is `(h_i,h_{i+m},h_{i+2m},\ldots)`.  This function returns the strip-decomposition of the Dyck path given by the height sequence 'h'.

	EXAMPLES::
	
		sage: [decompose(convert_to_height(x),2,3) for x in mTamari(2,3)]
		[[(1, 0, 1, 0, 1, 0, 0), (1, 0, 1, 0, 1, 0, 0)],
		 [(1, 0, 1, 0, 1, 0, 0), (1, 0, 1, 1, 0, 0, 0)],
		 [(1, 0, 1, 1, 0, 0, 0), (1, 0, 1, 1, 0, 0, 0)],
		 [(1, 0, 1, 0, 1, 0, 0), (1, 1, 0, 0, 1, 0, 0)],
		 [(1, 0, 1, 0, 1, 0, 0), (1, 1, 0, 1, 0, 0, 0)],
		 [(1, 0, 1, 1, 0, 0, 0), (1, 1, 0, 1, 0, 0, 0)],
		 [(1, 0, 1, 1, 0, 0, 0), (1, 1, 1, 0, 0, 0, 0)],
		 [(1, 1, 0, 0, 1, 0, 0), (1, 1, 0, 0, 1, 0, 0)],
		 [(1, 1, 0, 0, 1, 0, 0), (1, 1, 0, 1, 0, 0, 0)],
		 [(1, 1, 0, 1, 0, 0, 0), (1, 1, 0, 1, 0, 0, 0)],
		 [(1, 1, 0, 1, 0, 0, 0), (1, 1, 1, 0, 0, 0, 0)],
		 [(1, 1, 1, 0, 0, 0, 0), (1, 1, 1, 0, 0, 0, 0)]]

	REFERENCES::
	
		.. [KM15] Myrto Kallipoliti, Henri M\"uhle.
		   *The m-Cover Posets and Their Applications*
		   Adv. Appl. Math., *69*, pages 65-108, 2015.
	"""
	return [convert_from_height(q,1,n) for q in [tuple([h[i+j*m-1] for j in range(n)]) for i in range(1,m+1)]]
	
def decompose_all(m,n):
	"""
	Returns the list of strip-decompositions of all m-Dyck paths of parameter n.
	
	:params m: an integer.
	:params n: an integer.
	
	EXAMPLES::
	
		sage: decompose_all(2,3)
		[[(1, 0, 1, 0, 1, 0, 0), (1, 0, 1, 0, 1, 0, 0)],
		 [(1, 0, 1, 0, 1, 0, 0), (1, 0, 1, 1, 0, 0, 0)],
		 [(1, 0, 1, 1, 0, 0, 0), (1, 0, 1, 1, 0, 0, 0)],
		 [(1, 0, 1, 0, 1, 0, 0), (1, 1, 0, 0, 1, 0, 0)],
		 [(1, 0, 1, 0, 1, 0, 0), (1, 1, 0, 1, 0, 0, 0)],
		 [(1, 0, 1, 1, 0, 0, 0), (1, 1, 0, 1, 0, 0, 0)],
		 [(1, 0, 1, 1, 0, 0, 0), (1, 1, 1, 0, 0, 0, 0)],
		 [(1, 1, 0, 0, 1, 0, 0), (1, 1, 0, 0, 1, 0, 0)],
		 [(1, 1, 0, 0, 1, 0, 0), (1, 1, 0, 1, 0, 0, 0)],
		 [(1, 1, 0, 1, 0, 0, 0), (1, 1, 0, 1, 0, 0, 0)],
		 [(1, 1, 0, 1, 0, 0, 0), (1, 1, 1, 0, 0, 0, 0)],
		 [(1, 1, 1, 0, 0, 0, 0), (1, 1, 1, 0, 0, 0, 0)]]
	"""
	return [decompose(h,m,n) for h in convert_to_height_all(m,n)]

#######################
#####RENAMING##########
#######################
	
def rename_decomposed_paths(m,n):
	"""
	Renames the set of strip-decompositions of all m-Dyck paths of parameter n, where each occurring Dyck path is replaced by its index in a fixed linear extension of mTamari(1,n)
	
	:params m: an integer.
	:params n: an integer.
	"""
	D = decompose_all(m,n)
	L = mTamari(1,n).list()
	return [tuple([L.index(d[i]) for i in range(m)]) for d in D]
	
def rename_Tamari(T):
	"""
	Relabels the elements in the Tamari lattice by their index in a fixed linear extension.
	
	:params T: a Tamari lattice.
	"""
	return LatticePoset(T.relabel({T[i]:i for i in range(len(T))}))

#######################
#####BOUNCING##########
#######################

def bounce_pair(T,d,i,j):
	"""
	Applies the bounce map from Section~3.3 in [KM15]_ to the i-th and j-th entry in the m-tuple 'd' of Dyck paths.
	
	In an m-tuple `d` of Dyck paths, the bounce map `beta_{i,j}` replaces the i-th entry by the meet of d[i] and d[j], and the j-th entry by their join.  This function performs the action of `beta_{i,j}` on the given tuple `d`.  
	
	:params T: a Tamari lattice.
	:params d: an m-tuple of integers.
	:params i: an integer.
	:params j: an integer.

	EXAMPLES::
	
		sage: D = rename_decomposed_paths(3,3)
		sage: T = rename_Tamari(mTamari(1,3))
		sage: [bounce_pair(T,d,0,1) for d in D]
		[(0, 0, 0),
		 (0, 0, 1),
		 (0, 1, 1),
		 (1, 1, 1),
		 (0, 0, 2),
		 (0, 0, 3),
		 (0, 1, 3),
		 (1, 1, 3),
		 (1, 1, 4),
		 (0, 2, 2),
		 (0, 2, 3),
		 (0, 3, 3),
		 (0, 4, 3),
		 (0, 4, 4),
		 (1, 4, 4),
		 (2, 2, 2),
		 (2, 2, 3),
		 (2, 3, 3),
		 (3, 3, 3),
		 (3, 3, 4),
		 (3, 4, 4),
		 (4, 4, 4)]

	REFERENCES::
	
		.. [KM15] Myrto Kallipoliti, Henri M\"uhle.
		   *The m-Cover Posets and Their Applications*
		   Adv. Appl. Math., *69*, pages 65-108, 2015.
	"""
	p = d[i]
	q = d[j]
	L = len(d)*[0]
	for k in range(len(d)):
		if k==i:
			L[k] = T.meet(p,q)
		elif k==j:
			L[k] = T.join(p,q)
		else:
			L[k] = d[k]
	d = tuple(L)
	return d

def bounce(T,d):
	"""
	Applies the bounce map `beta_{i,j}` to the given m-tuple of Dyck paths for `1 \leq i < j \leq m` in that order.
	
	:params T: a Tamari lattice.
	:params d: an m-tuple of integers.

	EXAMPLES::
	
		sage: D = rename_decomposed_paths(3,3)
		sage: T = rename_Tamari(mTamari(1,3))
		sage: [bounce(T,d) for d in D]
		[(0, 0, 0),
		 (0, 0, 1),
		 (0, 1, 1),
		 (1, 1, 1),
		 (0, 0, 2),
		 (0, 0, 3),
		 (0, 0, 4),
		 (0, 1, 4),
		 (1, 1, 4),
		 (0, 2, 2),
		 (0, 2, 3),
		 (0, 3, 3),
		 (0, 3, 4),
		 (0, 4, 4),
		 (1, 4, 4),
		 (2, 2, 2),
		 (2, 2, 3),
		 (2, 3, 3),
		 (3, 3, 3),
		 (3, 3, 4),
		 (3, 4, 4),
		 (4, 4, 4)]
	"""
	for i in range(len(d)):
		for j in range(i+1,len(d)):
			d = bounce_pair(T,d,i,j)
	return d

def rename_and_bounce(m,n):
	"""
	Applies the bounce map to the m-Tamari lattice of parameter n.
	
	:params m: an integer.
	:params n: an integer.
	"""
	Tr = rename_Tamari(mTamari(1,n))
	return [bounce(Tr,d) for d in rename_decomposed_paths(m,n)]
	
#######################
#####CHECKING##########
#######################
	
def componentwise_subposet(T,L):
	"""
	Constructs the subposet of T^m induced by the elements in the list L.
	
	:param T: a Tamari lattice.
	:param L: a list of m-tuples of integers.

	EXAMPLES::
	
		sage: componentwise_subposet(mTamari(1,3),rename_and_bounce(3,3))
		Finite poset containing 22 elements
	"""
	Tr = rename_Tamari(T)
	P = []
	for l in L:
		C = []
		for k in L:
			if k==l:
				continue
			if all([Tr.le(l[i],k[i]) for i in range(len(l))]):
				C = C + [k]
		P = P + [C]
	D = dict([[L[i],P[i]] for i in range(len(L))])
	return Poset(D)

def check_bounce_conjecture(m1,m2,n1,n2):
	"""
	Checks whether Conjecture~3.28 of [KM15]_ holds for all `m1 \leq m \leq m2` and all `n1 \leq n \leq n2`.
	
	Conjecture~3.28 of [KM15]_ claims that the m-Tamari lattice of parameter n is isomorphic to the componentwise subposet of the m-fold direct product of the Tamari lattice of parameter n induced by the bounced strip-decomposition of all m-Dyck paths of parameter n.  This function checks whether this conjecture is true for m and n in the given range.
	
	:params m1: an integer.
	:params m2: an integer.
	:params n1: an integer.
	:params n2: an integer.
	
	EXAMPLES::
	
		sage: check_bounce_conjecture(1,2,3,5)
		n = 3 / m = 1 :
			Tamari lattice created! ( 0.0 mins )
			Decomposition done! ( 0.0 mins )
			Poset created! ( 0.0 mins )
			Success! ( 0.0 mins )
		n = 3 / m = 2 :
			Tamari lattice created! ( 0.0 mins )
			Decomposition done! ( 0.0 mins )
			Poset created! ( 0.0 mins )
			Success! ( 0.0 mins )
		-----------------------
		n = 4 / m = 1 :
			Tamari lattice created! ( 0.0 mins )
			Decomposition done! ( 0.0 mins )
			Poset created! ( 0.0 mins )
			Success! ( 0.0 mins )
		n = 4 / m = 2 :
			Tamari lattice created! ( 0.0 mins )
			Decomposition done! ( 0.0 mins )
			Poset created! ( 0.0 mins )
			Success! ( 0.0 mins )
		-----------------------
		n = 5 / m = 1 :
			Tamari lattice created! ( 0.0 mins )
			Decomposition done! ( 0.0 mins )
			Poset created! ( 0.0 mins )
			Success! ( 0.0 mins )
		n = 5 / m = 2 :
			Tamari lattice created! ( 0.0 mins )
			Decomposition done! ( 0.0 mins )
			Poset created! ( 0.01 mins )
			Success! ( 0.0 mins )
		-----------------------

	REFERENCES::
		.. [KM15] Myrto Kallipoliti, Henri M\"uhle.
		   *The m-Cover Posets and Their Applications*
		   Adv. Appl. Math., *69*, pages 65-108, 2015.
	"""
	import time
	for i in range(n1,n2+1):
		for j in range(m1,m2+1):
			a0 = time.time()
			print 'n =',i,'/ m =',j,':'
			T = mTamari(1,i)
			a1 = time.time()
			print '\t Tamari lattice created! (',round((a1-a0)/60,2),'mins )'
			L = rename_and_bounce(j,i)
			a2 = time.time()
			print '\t Decomposition done! (',round((a2-a1)/60,2),'mins )'
			P = componentwise_subposet(T,L)
			a3 = time.time()
			print '\t Poset created! (',round((a3-a2)/60,2),'mins )'
			if P.is_isomorphic(mTamari(j,i)):
				a4 = time.time()
				print '\t Success! (',round((a4-a3)/60,2),'mins )'
			else:
				a4 = time.time()
				print '\t Failure! (',round((a4-a3)/60,2),'mins )'
		print '-----------------------'