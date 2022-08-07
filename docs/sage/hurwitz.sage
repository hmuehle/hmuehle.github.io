def deduce_all_words(word,rules,max_counter=100):
	# input: word .. string; rules .. a list of rules; max_counter .. integer
	#  each letter in 'word' is considered as a generator
	#  each rule in 'rules' is a list of two-letter strings 
	# output: a list containing all factorizations in the same Hurwitz orbit of 'word', and the corresponding Hurwitz graph
	n = len(word)
	new_words = [word]
	used_words = []
	new_cycles = dict()
	last_change = 0
	hurwitz_graph = Graph()
	counter = 0
	suggestions = dict()

	while true:
		##################
		#create new words#
		##################
		act_word = new_words.pop(0)
		change_count = 0            #logs whether a new word was found
		for i in range(n-1):
			pr = act_word[i:i+2]
			for r in rules:
				if pr in r:
					change_count += 1
					for repl in r:
						substitute = act_word[0:i] + repl + act_word[i+2:n]
						if act_word==substitute:
							continue
						if not hurwitz_graph.has_edge(act_word,substitute):
							hurwitz_graph.add_edge(act_word,substitute)
						if substitute not in new_words and substitute not in used_words and substitute!=act_word:
							new_words.append(substitute)
					#there can be only one rule applicable to the 2-letter subword starting at position i
					break
		last_change += 1
		if change_count == n-1:
			#if we have used all possible Hurwitz-moves, we are done with this word
			used_words.append(act_word)
			last_change = 0
		else:
			new_words.append(act_word)
		###################
		#detect new cycles#
		###################
		# i.e. for each word in 'new_words' we check for length-2 subwords that do not appear in rules
		# if we find such a word and its reverse already appears in 'rules', then we know the length of the Hurwitz orbit
		for test_word in new_words:
			for i in range(n-1):
				pr = test_word[i:i+2]
				sug_length = 0
				new = true
				for r in rules:
					if pr in r:
						new = false
						break
				if new and pr not in new_cycles:
					if pr[0:1]==pr[1:2]:
						if tuple([pr]) in suggestions.keys():
							suggestions[tuple([pr])].append(test_word)
						else:
							suggestions[tuple([pr])] = [test_word]
						sug_length = 1
					for r in rules:
						rp = pr[1:2] + pr[0:1]
						if rp in r:
							sug_length = len(r)
					new_cycles[pr] = sug_length
		#####################################
		#determine if you're stuck in a loop#
		#####################################
		# this basically happens as soon as you've cycled once through 'new_words' and haven't completed a word
		# other indication is that you've run this loop for as many times as 'max_counter' suggests
		counter += 1
		if last_change>len(new_words) or counter==max_counter:
			####################################
			#check if you can suggest relations#
			####################################
			# if there are two words 'yabz' and 'ycdz', then we know that ab=cd
			test_words = new_words + used_words
			for i in range(len(test_words)):
				word_1 = test_words[i]
				for j in range(i+1,len(test_words)):
					word_2 = test_words[j]
					for k in range(n-1):
						x1 = word_1[:k]
						y1 = word_1[k+2:]
						x2 = word_2[:k]
						y2 = word_2[k+2:]
						if x1==x2 and y1==y2:
							r1 = word_1[k:k+2]
							r2 = word_2[k:k+2]
							new_rule = true
							for r in rules:
								if r1 in r and r2 in r:
									new_rule = false
									break
							if new_rule:
								tt_1 = tuple([r1,r2])
								tt_2 = tuple([r2,r1])
								if tt_1 not in suggestions.keys() and tt_2 not in suggestions.keys():
									suggestions[tt_1] = [[word_1,word_2]]
								else:
									if tt_1 in suggestions.keys():
										suggestions[tt_1].append([word_1,word_2])
									if tt_2 in suggestions.keys():
										suggestions[tt_2].append([word_1,word_2])
			if counter==max_counter:
				message = 'Traveled a long way. '
			else:
				message = 'Trapped in a loop! '
			print message, 'Try adding new rules involving the new cycles:',new_cycles
			print 'New words: (',len(new_words),') ',new_words
			print 'Used words: (',len(used_words),')',used_words
			print 'Suggestions:'
			for sug in suggestions.keys():
				print '\t',sug,suggestions[sug]
			#############################
			#check for unnecessary rules#
			#############################
			# if there is a rule containing a word 'ab' that does not appear anywhere in 'new_words' or 'used_words', then we can omit it
			bad_rules = []
			for r in rules:
				good_rule = false
				for w in new_words+used_words:
					if w.find(r[0])>-1:
						good_rule = true
						break
				if not good_rule:
					bad_rules.append(r)
			print 'Unnecessary Rules:',bad_rules
			break
		######################
		#check if you're done#
		######################
		if len(new_words)==0:
			break
	return [used_words,hurwitz_graph]

def create_factorization_poset(words,rules,max_counter=100):
	# input: words .. a complete Hurwitz orbit, i.e. a list of strings; rules .. a list of rules; max_counter .. integer
	# output: the factorization poset whose maximal chains correspond to 'words, and a dictionary of equivalent factorizations
	######################
	#collect all elements#
	######################
	top = words[0]
	n = len(top)
	equivalent_words = {}
	equivalent_words[''] = ['']
	equivalent_words[top] = words
	for r in rules:
		equivalent_words[r[0]] = r
		for x in r:
			l = x[0:1]
			if l not in equivalent_words.keys():
				equivalent_words[l] = [l]
	for k in range(3,n):
		for w in words:
			k_prefix = w[0:k]
			if len(equivalent_words.keys())==0:
				equivalent_words[k_prefix] = deduce_all_words(k_prefix,rules,max_counter)[0]
			else:
				known = false
				for x in equivalent_words.keys():
					if k_prefix in equivalent_words[x]:
						known = true
						break
				if not known:
					equivalent_words[k_prefix] = deduce_all_words(k_prefix,rules,max_counter)[0]
	#######################################################
	#creating the poset diagram of the factorization poset#
	#######################################################
	H = DiGraph()
	for x1 in equivalent_words.keys():
		for x2 in equivalent_words.keys():
			if len(x1)==len(x2)-1:
				for w in equivalent_words[x2]:
					if x1 == w[:len(x1)]:
						H.add_edge(x1,x2)
						break
	return [Poset(H),equivalent_words]

def check_common_cover(P,order):
	# input: P .. a poset; order .. a total order of the atoms of 'P'
	# output: true if and only if 'P' satisfies the common cover property with respect to 'order'
	#  in other words: it checks if for each non-minimal element x in 'order' there exists some element y preceding x such that x and y have a common upper cover
	for i in range(1,len(order)):
		x = order[i]
		val = false
		for j in range(i):
			y = order[j]
			if len(set(P.upper_covers(x)) & set(P.upper_covers(y)))>0:
				val = true
				break
		if not val:
			print 'Counterexample!', x
			return false
	return true

def extract_cycles(facts):
	# input: equivalent_words .. dictionary of equivalent factorizations of prefixes of some element
	# output: a list of Hurwitz orbits of prefixes of length 2
	cycles = []
	for f in facts.keys():
		if len(f)==2:
			cycles.append(facts[f])
	return cycles

def all_tuples(val):
	# input: val .. list of integers
	# output: all tuples of length 'len(val)' with entries in {0,1,...,'val[i]'} for all i
	ranges = []
	for i in val:
		ranges.append(range(i))
	return cartesian_product(ranges)

def check_for_compatible_order(cycles,Q):
	# input: cycles .. list of Hurwitz orbits of rank-2 elements in a factorization poset; 'Q' .. the corresponding factorization poset
	# output: lists all compatible orders
	cycle_sizes = []
	for l in cycles:
		cycle_sizes.append(len(l))
	print cycle_sizes
	k = all_tuples(cycle_sizes)
	print 'Computed all tuples!'
	c = 0
	cnt = 0
	tot = prod(cycle_sizes)
	fails = []
	for t in k:
		cnt += 1
		print '\rProgress:',cnt,'/',tot,
		if true:
			D = DiGraph()
			for i in range(len(tuple(t))):
				if len(cycles[i])==1:
					continue
				for j in range(len(cycles[i])):
					if j==t[i]:
						D.add_edge(cycles[i][j][0:1],cycles[i][j][1:2])
					else:
						D.add_edge(cycles[i][j][1:2],cycles[i][j][0:1])
			if D.is_directed_acyclic():
				P = Poset(D)
				if not P.is_chain():
					print t,'Something is wrong here!!'
				else:
					c += 1
					q = P.maximal_chains()[0]
					val = check_common_cover(Q,q)
					print t,q,val
					if not val:
						fails.append(q)
	print '\nDone.  ' + str(c) + ' results.'
	return fails