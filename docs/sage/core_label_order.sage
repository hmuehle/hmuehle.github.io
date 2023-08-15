def jsd_label(L,u,v):
	r"""
	Returns the smallest element of L whose join with u is v.
	
	:param L: a finite lattice
	:param u: an element of L
	:param v: an element of L
	
	EXAMPLES::

		sage: L = LatticePoset(DiGraph({0:[1,2],1:[3],2:[4],3:[4]}))
		sage: jsd_label(L,3,4)
		2
	"""
	return L.meet([x for x in L if L.join(u,x)==v])

def jsd_labeling(L):
	r"""
	Returns the canonical edge-labeling of a join-semidistriibutive lattice.  This is the map \eta used in Proposition 3.4 of [Bar19]
	
	:param L: a finite join-semidistributive lattice
	
	EXAMPLES::

		sage: L = LatticePoset(DiGraph({0:[1,2],1:[3],2:[4],3:[4]}))
		sage: labeling = jsd_labeling(L)
		sage: for cov in L.cover_relations():
		....:     cov,labeling(cov[0],cov[1])
		....:
		([0, 1], 1)
		([0, 2], 2)
		([1, 3], 3)
		([3, 4], 2)
		([2, 4], 1)
		
	REFERENCES::
		.. [Bar19] Emily Barnard
		   *The Canonical Join Complex*
		   The Electronic Journal of Combinatorics (26), 2019.
	"""
	assert(L.is_join_semidistributive())
	
	return lambda u,v:jsd_label(L,u,v)

def trim_label(L,X,u,v):
	r"""
	Returns the index i of the first element x_i of X such that u\vee(x_i\wedge v)=v. 
	
	:param L: a finite lattice
	:param X: a maximal chain of L
	:param u: an element of L
	:param v: an element of L
	
	EXAMPLES::

		sage: L = LatticePoset(DiGraph({0:[1,2],1:[3],2:[4],3:[4]}))
		sage: X = [0,1,3,4]
		sage: trim_label(L,X,3,4)
		3
	"""
	val = -1
	for i in range(len(X)):
		if L.join(u,L.meet(X[i],v))==v:
			val = i
			break
	return val

def trim_labeling(L):
	r"""
	Returns the interpolating labeling of a trim lattice.  This is the labeling \gamma_3 appearing in Proposition 2 of [Tho06].

	:param L: a finite trim lattice

	EXAMPLES::

		sage: L = LatticePoset(DiGraph({0:[1,2],1:[3],2:[4],3:[4]}))
		sage: labeling = trim_labeling(L)
		sage: for cov in L.cover_relations():
		....:     cov,labeling(cov[0],cov[1])
		....:
		([0, 1], 1)
		([0, 2], 3)
		([1, 3], 2)
		([3, 4], 3)
		([2, 4], 1)

	REFERENCES::
		.. [Tho06] Hugh Thomas
		   *An Analogue of Distributivity for Ungraded Lattices*
		   Order (23), 2006.
	"""
	assert(L.is_trim())
	
	X = [C for C in L.maximal_chains() if len(C)==L.rank()+1][0]
	return lambda u,v:trim_label(L,X,u,v)

def core_label_sets(L,labeling,dual=false):
	r"""
	Returns the core label sets of L with respect to the given edge-labeling.
	
	The *core label set* of am element x in a finite lattice is the set of edge labels appearing in the interval [x',x], where x'=0\wedge x_1\wedge x_2\wedge\cdots\wedge x_k, where 0 is the bottom element of L and \{x_1, x_2, ..., x_k\} is the set of lower covers of x.  See Section~3.1 in [Müh22].

	:param L: a finite lattice
	:param labeling: an edge-labeling of L
	:param dual: a boolean parameter toggling whether core label sets are determined through the meet of the lower covers of x or through the join of the upper covers of x

	EXAMPLES::

		sage: L = LatticePoset(DiGraph({0:[1,2],1:[3],2:[4],3:[4]}))
		sage: core_label_sets(L,jsd_labeling(L))
		{0: {}, 1: {1}, 3: {3}, 2: {2}, 4: {1, 2, 3}}

	REFERENCES::
		.. [Müh22] Henri Mühle
		   *Meet-Distributive Lattices have the Intersection Property*
		   Mathematica Bohemica (148), 2022.
	"""
	dl_sets = dict()
	for x in L:
		xx = -1
		QQ = -1
		if dual:
			xx = L.join(L.upper_covers(x))
			QQ = L.interval(x,xx)
		else:
			xx = L.meet(L.lower_covers(x))
			QQ = L.interval(xx,x)
		labs = []
		for cov in L.cover_relations():
			if cov[0] in QQ and cov[1] in QQ:
				labs.append(labeling(cov[0],cov[1]))
		dl_sets[x] = Set(labs)
	return dl_sets

def core_label_order(L,labeling,dual=false):
	r"""
	Returns the core label order of L.
	
	The *core label order* is the alternative order on L, determined through containment of its core label sets.  See Section~3.1 in [Müh22].
	:param L: a finite lattice
	:param labeling: an edge labeling of L
	:param dual: a boolean parameter toggling whether core label sets are determined through the meet of the lower covers of x or through the join of the upper covers of x

	EXAMPLES::

		sage: L = LatticePoset(DiGraph({0:[1,2],1:[3],2:[4],3:[4]}))
		sage: C = core_label_order(L,jsd_labeling(L)); C
		Finite poset containing 5 elements
		sage: C.cover_relations()
		[[{}, {3}],
		 [{}, {2}],
		 [{}, {1}],
		 [{3}, {1, 2, 3}],
		 [{2}, {1, 2, 3}],
		 [{1}, {1, 2, 3}]]
	
	REFERENCES::
		.. [Müh22] Henri Mühle
		   *Meet-Distributive Lattices have the Intersection Property*
		   Mathematica Bohemica (148), 2022.
	"""
	G = DiGraph()
	C = core_label_sets(L,labeling,dual)
	for c in C.keys():
		G.add_vertex(C[c])
	for c1 in G.vertices():
		for c2 in G.vertices():
			if c1==c2:
				continue
			if c1.issubset(c2):
				G.add_edge(c1,c2)
	return Poset(G)

def is_core_labeling(L,labeling):
	r"""
	Returns whether the given edge-labeling is a core labeling.
	
	A *core labeling* is an edge-labeling of L such that no two elements of L have the same core label set.  See Section~3.1 in [Müh22].

	:param L: a finite lattice
	:param labeling: an edge-labeling

	EXAMPLES::

		sage: L = LatticePoset(DiGraph({0:[1,2],1:[3],2:[4],3:[4]}))
		sage: is_core_labeling(L,jsd_labeling(L))
		True
		sage: labeling = lambda u,v:1
		sage: is_core_labeling(L,labeling)
		False

	REFERENCES::
		.. [Müh22] Henri Mühle
		   *Meet-Distributive Lattices have the Intersection Property*
		   Mathematica Bohemica (148), 2022.
	"""
	cl_set = core_label_sets(L,labeling)
	for i in cl_set.keys():
		C1 = cl_set[i]
		for j in cl_set.keys():
			if i==j:
				continue
			C2 = cl_set[j]
			if C1.issubset(C2) and C2.issubset(C1):
				return false
	return true
	
def canonical_join_representation(L,u):
	r"""
	Returns the canonical join representation of u.
	
	The *canonical join representation* of u is a minimal antichain of join-irreducible elements whose join is u.  See Section I.3 in [Fre95].

	:param L: a finite lattice
	:param u: an element of L

	EXAMPLES::

		sage: L = LatticePoset(DiGraph({0:[1,2],1:[3],2:[4],3:[4]}))
		sage: canonical_join_representation(L,4)
		[2, 1]
		
	REFERENCES::
		.. [Fre95] Ralph Freese, Jaroslav Jezek, James B. Nation
		   *Free Lattices*
		   American Mathematical Society, 1995.
	"""
	return [jsd_label(L,v,u) for v in L.lower_covers(u)]
	
def kappa_map(L,u):
	r"""
	Returns the image of the kappa map applied to u.
	
	The *kappa map* assigns to a given join-irreducible element j of a lattice L the largest element in L (if it exists) which is below the unique lower cover of j, but not below j itself. This is a bijection between join- and meet-irreducible elements when L is semidistributive.  See Theorem~2.54 in [Fre95].

	:param L: a finite lattice
	:param u: a join-irreducible element of L

	EXAMPLES::

		sage: L = LatticePoset(DiGraph({0:[1,2],1:[3],2:[4],3:[4]}))
		sage: [(j,kappa_map(L,j)) for j in L.join_irreducibles()]
		[(1, 2), (3, 1), (2, 3)]
		sage: M = LatticePoset(DiGraph({0:[1,2,3],1:[4],2:[4],3:[4]}))
		sage: [(j,kappa_map(M,j)) for j in M.join_irreducibles()]
		[(1, False), (2, False), (3, False)]

	REFERENCES::
		.. [Fre95] Ralph Freese, Jaroslav Jezek, James B. Nation
		   *Free Lattices*
		   American Mathematical Society, 1995.
	"""
	assert(u in L.join_irreducibles())

	v = L.lower_covers(u)[0]
	candidates = [q for q in L if L.le(v,q) and not L.le(u,q)]
	C = L.subposet(candidates)
	if C.has_top():
		return C.top()
	else:
		return false

def extended_kappa_map(L,u):
	r"""
	Returns the image of the extended kappa map applied to u.
	
	The *extended kappa map* assigns to a given element u meet of the images of the kappa map applied to the canonical join representation of u.  See Definition 1.1.3 in [Bar21].

	:param L: a finite lattice
	:param u: an element of L

	EXAMPLES::

		sage: L = LatticePoset(DiGraph({0:[1,2],1:[3],2:[4],3:[4]}))
		sage: [(u,extended_kappa_map(L,u)) for u in L]
		[(0, 4), (1, 2), (3, 1), (2, 3), (4, 0)]
	
	REFERENCES::
		.. [Bar21] Emily Barnard, Gordana Todorov, Shijie Zhu
		   *Dynamical Combinatorics and Torsion Classes*
		   Journal of Pure and Applied Algebra (225), 2021.
	"""
	assert(L.is_join_semidistributive)

	return L.meet([kappa_map(L,j) for j in canonical_join_representation(L,u)])

def kappa_leq(L,u,v):
	r"""
	Returns whether u is less than or equal to v with respect to the kappa order.

	:param L: a finite lattice
	:param u: an element of L
	:param v: an element of L

	EXAMPLES::

		sage: L = LatticePoset(DiGraph({0:[1,2],1:[3],2:[4],3:[4]}))
		sage: kappa_leq(L,2,4)
		True
		sage: kappa_leq(L,2,1)
		False
	"""
	assert(L.is_join_semidistributive())

	return L.le(u,v) and L.le(extended_kappa_map(L,v),extended_kappa_map(L,u))

def kappa_order(L):
	r"""
	Returns the kappa order of a semidistributive lattice.
	
	The *kappa order* is the alternative order on a semidistributive lattice L which is determined by comparing u with v in L and also comparing the images of v and u under the extended kappa map.  See Definition 4.21 in [Eno23].

	:param L: a finite semidistributive lattice

	EXAMPLES::

		sage: L = LatticePoset(DiGraph({0:[1,2],1:[3],2:[4],3:[4]}))
		sage: K = kappa_order(L); K
		Finite poset containing 5 elements
		sage: K.cover_relations()
		[[0, 1], [0, 2], [0, 3], [1, 4], [2, 4], [3, 4]]
		
	REFERENCES::
		.. [Eno23] Haruhisa Enomoto
		   *From the Lattice of Torsion Classes to the Posets of Wide Subcategories and ICE-closed Subcategories*
		   Algebras and Representation Theory, 2023.
	"""
	assert(L.is_semidistributive())
	
	order = lambda u,v:kappa_leq(L,u,v)
	return Poset([list(L),order])
