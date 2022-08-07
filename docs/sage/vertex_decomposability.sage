##########################
###SIMPLICIAL COMPLEXES###
##########################
	
def diagonal_complex(L):
	"""
	Constructs a simplicial complex whose facets are the elements of the given list.
	
	:param L: a list.
	
	EXAMPLES::
	
		sage: diagonal_complex(range(3))
		Simplicial complex with vertex set (0, 1, 2) and facets {(2,), (0,), (1,)}
	"""
	return SimplicialComplex([[x] for x in L])
    
def simplex(n):
	"""
	Constructs a simplex with n vertices.
	
	A *simplex* is a simplicial complex with a unique facet.
	
	:param n: an integer.
	
	EXAMPLES::
	
		sage: simplex(3)
		Simplicial complex with vertex set (0, 1, 2) and facets {(0, 1, 2)}
	"""
	return SimplicialComplex([range(n)])
    
def is_simplex(C):
	"""
	Checks whether the given simplicial complex is a simplex, i.e. if it has a unique facet.
	
	:param C: a simplicial complex.
	
	EXAMPLES::
	
		sage: is_simplex(SimplicialComplex([[1,2],[1,3],[2,3]]))
		False
		sage: is_simplex(simplex(3))
		True

	"""
	return len(C.facets())==1
	
############################
###VERTEX DECOMPOSABILITY###
############################

def delete_vertex(C,v):
	"""
	The deletion of a vertex in a given simplicial complex.
	
	The *deletion* of a vertex v is the simplicial complex formed by all faces that do not contain v.
	
	:param C: a simplicial complex.
	:param v: a vertex.
	
	EXAMPLES::
	
		sage: C = SimplicialComplex([[1,2],[1,3],[2,3]])
		sage: delete_vertex(C,1)
		Simplicial complex with vertex set (2, 3) and facets {(2, 3)}
	"""
	if v not in C.vertices():
		return C
	fcs = C.faces()
	new_fcs = dict()
	for i in fcs.keys():
		L = []
		for X in fcs[i]:
			if v not in X:
				L.append(X)
		new_fcs[i] = L
	
	N = max(new_fcs.keys())
	facets = new_fcs[N]
	for i in range(-1,N):
		for X in new_fcs[i]:
			is_max = true
			for Y in new_fcs[i+1]:
				if Set(X).issubset(Set(Y)):
					is_max = false
					break
			if is_max:
				facets.append(X)
	return SimplicialComplex(facets)                                    

def is_shedding_vertex(C,v):
	"""
	Checks whether a vertex is a shedding vertex of a simplicial complex.
	
	A vertex `v` of a simplicial complex `\Delta` is a *shedding vertex* if the following three conditions hold:
	
	1. The link of `v` in `\Delta` is vertex-decomposable.
	
	2. The deletion of `v` in `\Delta` is vertex-decomposable.
	
	3. The link of `v` in `\Delta` and the deletion of `v` in `\Delta` do not share a common facet.
	
	:param C: a simplicial complex.
	:param v: a vertex.
	
	EXAMPLES::
	
		sage: C = SimplicialComplex([[1,2],[1,3],[2,3]])
		sage: is_shedding_vertex(C,1)
		True
		
		sage: D = SimplicialComplex([[1,2],[2,3],[3,4]])
		sage: is_shedding_vertex(D,2)
		False
	"""
	lk = C.link([v])
	rm = delete_vertex(C,v)
	if Set(lk.facets()).intersection(Set(rm.facets())).cardinality()>0:
		return false
	if not is_vertex_decomposable(lk):
		return false
	if not is_vertex_decomposable(rm):
		return false
	return true
	
def shedding_vertices(C):
	"""
	Computes the set of shedding vertices of a simplicial complex.
	
	:param C: a simplicial complex.
	
	EXAMPLES::
	
		sage: C = SimplicialComplex([[1,2],[1,3],[2,3]])
		sage: shedding_vertices(C)
		[1, 2, 3]
		sage: D = SimplicialComplex([[1,2],[2,3],[3,4]])
		sage: shedding_vertices(D)
		[1, 4]
	"""
	
	if is_simplex(C):
		return list(C.vertices())
	sheds = []
	for v in C.vertices():
		lk = C.link([v])
		rm = delete_vertex(C,v)
		if Set(lk.facets()).intersection(Set(rm.facets())).cardinality()>0:
			continue
		if not is_vertex_decomposable(lk):
			continue
		if not is_vertex_decomposable(rm):
			continue
		sheds.append(v)
	return sheds

def is_vertex_decomposable(C):
	"""
	Checks whether a simplicial complex is vertex-decomposable.
	
	A simplicial complex is *vertex-decomposable*, if it is either empty, or a simplex, or if it has a shedding vertex.
	
	:param C: a simplicial complex.
	
	EXAMPLES::
		sage: C = SimplicialComplex([[1,2],[1,3],[2,3]])
		sage: is_vertex_decomposable(C)
		True
		sage: D = SimplicialComplex([[1,2],[2,3],[3,4]])
		sage: is_vertex_decomposable(D)
		True
		sage: E = SimplicialComplex([[1,2],[3,4]])
		sage: is_vertex_decomposable(E)
		False
	"""
	if is_simplex(C):
		return true
	else:
		return len(shedding_vertices(C))>0