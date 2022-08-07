######################
###SPERNER PROPERTY###
######################

def rank_vector(P):
	"""
	Returns the rank vector of a ranked poset.
	
	The *rank vector* of a ranked poset is the list of rank sizes.
	
	:param P: a ranked poset.
	
	EXAMPLES::
	
		sage: P = posets.BooleanLattice(4)
		sage: rank_vector(P)
		[1, 4, 6, 4, 1]
	"""
	assert(P.is_ranked())
	rv = (P.rank()+1)*[0]
	for p in P:
		rv[P.rank(p)] += 1
	return rv
	
def rank_elements(P):
	r"""
	Returns a list of the ranks of P.
	
	This function computes a list, whose `i`-th contains a list of all elements of rank `i` in the poset `P`.
	
	:param P: a ranked poset.
	
	EXAMPLES::
	
		sage: P = Poset({0:[3],1:[3],2:[4,5]})
		sage: rank_elements(P)
		[[2, 1, 0], [4, 5, 3]]
	"""
	assert(P.is_ranked())
	ranks = []
	for i in range(P.rank()+1):
		ranks.append([p for p in P if P.rank(p)==i])
	return ranks
	
def remove_rank(P,i):
	r"""
	Removes the the i-th rank from a poset.
	
	This function creates the subposet of `P` that consists of all elements, except those that have rank `i`.
	
	:param P: a ranked poset.
	:param i: an integer.
	
	EXAMPLES::
	
		sage: P = posets.BooleanLattice()
		sage: Q = remove_rank(P,1); Q
		Finite poset containing 12 elements
	"""
	assert(P.is_ranked())
	ranks = rank_elements(P)
	new_covs = dict()
	if i==0:
		for j in range(1,len(ranks)):
			for p in ranks[j]:
				new_covs[p] = P.upper_covers(p)
	elif i==P.rank():
		for j in range(len(ranks)-2):
			for p in ranks[j]:
				new_covs[p] = P.upper_covers(p)
		for p in ranks[i-1]:
			new_covs[p] = []
	else:
		for j in range(i-1):
			for p in ranks[j]:
				new_covs[p] = P.upper_covers(p)
		for j in range(i+1,len(ranks)-1):
			for p in ranks[j]:
				new_covs[p] = P.upper_covers(p)
		for p in ranks[i-1]:
			new_covs[p] = []
			for q in ranks[i+1]:
				if P.le(p,q):
					new_covs[p] += [q]
	return Poset(new_covs)
	
def is_sperner(P):
	"""
	Checks whether the given poset is Sperner.
	
	A ranked poset is *Sperner* if its width equals the size of its biggest rank.
	
	:param P: a ranked poset.
	
	EXAMPLES::
	
		sage: P = posets.BooleanLattice(4)
		sage: is_sperner(P)
		True
		
		sage: P = Poset({0:[3],1:[3],2:[4,5]})
		sage: is_sperner(P)
		False
	
	"""
	if P.cardinality()==0:
		return true
	return P.width()==max(rank_vector(P))

def is_strongly_sperner(P):
	r"""
	Checks whether the given poset is strongly Sperner.
	
	A ranked poset is *`k`-Sperner* if the sum of the sizes of the `k` largest ranks equals the sum of the sizes of the `k` largest antichains.  A ranked poset is *strongly Sperner* if it is `k`-Sperner for all `k`.  This property can be established by repeatedly removing the largest rank from the poset, and checking whether the remaining poset is 1-Sperner, until we obtain the empty poset.  See also Proposition 4.4 in [Mue17]_.
	
	:param P: a ranked poset.
	
	EXAMPLES::
	
		sage: P = Poset({0:[3],1:[4],2:[5,6],3:[7],4:[7],5:[8],6:[9]})
		sage: is_sperner(P)
		True
		sage: is_strongly_sperner(P)
		False

		sage: P = posets.BooleanLattice(4)
		sage: is_strongly_sperner(P)
		True
	
	REFERENCES::
	
		.. [Mue17] Henri Mühle.
		   *Symmetric Decompositions and the Strong Sperner Property for Noncrossing Partition Lattices*
		   J. Alg. Combin., *45*, pages 745-775. 2017.
	"""
	assert(P.is_ranked())
	val = []
	val += [is_sperner(P)]
	n = P.rank()
	Q = P.subposet(P.list())
	for k in range(1,n+1):
		rv = rank_vector(Q)
		l = 0
		# find maximal rank
		for i in range(len(rv)):
			if rv[i]>rv[l]:
				l = i
		# remove largest rank from poset
		Q = remove_rank(Q,l)
		val += [is_sperner(Q)]
	return prod(val)==1
	
###############################
####LYM PROPERTY###############
##(Lubell-Yamamoto-Meshalkin)##
###############################

def lym_value(P,A):
	r"""
	Computes the LYM-value of the given poset with respect to the given antichain.
	
	Let `P` be a ranked poset, and suppose that `N_i` is the number of elements of `P` of rank `i`.  Suppose that the given antichain `A` consists of `p_i` members of rank `i`.  The *LYM-value* (named after Lubell-Yamamoto-Meshalkin) is the sum of p_i/N_i.
	
	:param P: a ranked poset.
	:param A: an antichain of P.

	EXAMPLES::
	
		sage: P = Poset({0:[3],1:[4],2:[5,6],3:[7],4:[7],5:[8],6:[9]})
		sage: A = [5, 6, 1, 3] 
		sage: lym_value(P,A)
		13/12

		sage: P = posets.BooleanLattice(4)
		sage: A = [1, 2, 4]
		sage: lym_value(P,A)
		3/4
	"""
	
	assert(P.is_ranked())
	rv = rank_vector(P)
	v = (P.rank()+1)*[0]
	for a in A:
		v[P.rank(a)] += 1
	return sum([v[i]/rv[i] for i in range(P.rank()+1)])

def check_lym_poset(P):
	"""
	Checks whether the poset has a LYM-value not exceeding 1 for all its antichains.
	
	This function checks whether the LYM-value of the given poset with respect to each of its antichains is not exceeding 1.  We can discard antichains of size 2 or less, since their LYM-value can never be bigger than 1.  A poset with this property is strongly Sperner.
	
	:param P: a ranked poset.
	
	EXAMPLES:
	
		sage: P = Poset({0:[3],1:[4],2:[5,6],3:[7],4:[7],5:[8],6:[9]})
		sage: check_lym_poset(P)
		False
		
		sage: P = posets.BooleanLattice(4)
		sage: check_lym_poset(P)
		True
	"""
	
	assert(P.is_ranked())
	it = P.antichains_iterator()
	while True:
		try:
			A = it.next()
			if len(A)<=2:
				continue
			if lym_value(P,A)>1:
				return false
		except StopIteration:
			break
	return true

####################################
####NORMALIZED MATCHING PROPERTY####
####################################
	
def shade(P,A):
	"""
	Computes the shade of an antichain in a poset.
	
	Let `A` be a subset of the elements of `P` having rank `k`.  The *shade* of `A` consists of all elements of the poset that cover at least one member of this antichain.
	
	:param P: a ranked poset.
	:param A: an antichain of P.
	
	EXAMPLES::
	
		sage: P = Poset({0:[3],1:[4],2:[5,6],3:[7],4:[7],5:[8],6:[9]})
		sage: A = [5, 6, 3] 
		sage: shade(P,A)
		[8, 9, 7]

		sage: P = posets.BooleanLattice(4)
		sage: A = [1, 2, 4]
		sage: shade(P, A)
		[9, 3, 5, 10, 6, 12]
	"""
	assert(P.is_ranked())
	S = []
	for a in A:
		for p in P.upper_covers(a):
			if not p in S:
				S += [p]
	return S

def check_normalized_matching(P,A):
	r"""
	Checks whether the given poset has the normalized matching property with respect to the given antichain.
	
	Let `N_k` denote the number of elements of `P` of rank `k`.  A ranked poset has the *normalized matching property* if for every set `A` of elements of rank `k` the size of the shade of `A` divided by `N_{k+1}` is at least as big as the size of `A` divided by `N_k`.  
	
	:param P: a ranked poset.
	:param A: a list of elements of same rank.
	
	EXAMPLES::
	
		sage: P = Poset({0:[3],1:[4],2:[5,6],3:[7],4:[7],5:[8],6:[9]})
		sage: A = [1]
		sage: check_normalized_matching(P,A)
		False
	"""
	assert(P.is_ranked())
	if len(A)==0:
		return true
	k = P.rank(A[0])
	N = rank_vector(P)
	return len(shade(P,A))/N[k+1] >= len(A)/N[k]

def check_normalized_matching_all(P):
	"""
	Checks whether the normalized matching property holds for all antichains of fixed rank.
	
	This function checks whether for all k and for all subsets of elements of rank k the given poset has the normalized matching property.
	
	:param P: a ranked poset.

	EXAMPLES::
	
		sage: P = Poset({0:[3],1:[4],2:[5,6],3:[7],4:[7],5:[8],6:[9]})
		sage: check_normalized_matching_all(P)
		False
		
		sage: P = posets.BooleanLattice(4)
		sage: check_normalized_matching_all(P)
		True
	"""
	assert(P.is_ranked())
	ranks = rank_elements(P)
	#iterate through subsets of each rank, and check for the normalized matching property
	for k in range(P.rank()):
		rank_subsets = Subsets(ranks[k])
		for X in rank_subsets:
			if not check_normalized_matching(P,X):
				return false
	return true