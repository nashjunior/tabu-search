import argparse, sys, platform, os, operator, collections,random,sys
import numpy as np

path = os.path.dirname(os.path.abspath(__file__))
type_file = None
name = None
comment = None
dimension = None
edge_weight_type = None
edge_weight_format = None
display_data_type = None
distances = None

def arguments():
	parser = argparse.ArgumentParser(prog="Tabu search algorithm",usage='%(prog)s [options]')
	parser.add_argument('--file','-f',type=str,required=True,help='Arquivo a ser inserido')
	parser.add_argument('--maxiter','-mi',type=int,required=True,help='Quantidade de interacoes')
	parser.add_argument('--maxtabu','-mt',type=int,required=True,help='Quantidade de repeticoes proibidas')
	parser.add_argument('--maxcandidate','-mc',type=int,required=True,help='Quantidade maxima de candidados')
	args = parser.parse_args()
	global path
	if platform.system() == "Windows":
		path = path +r'\\'+args.file
	elif platform.system() == "Linux":
		path = path  + '/' + args.file
	decode_file(args.maxiter,args.maxtabu,args.maxcandidate)

def decode_file(maxiter,maxTabu,maxcandidate):
	global path, name, comment, dimension, edge_weight_type, edge_weight_format,distances
	with open(path,'r') as file:
		name = file.readline().strip().split()[1]
		type_file = file.readline().strip().split()[1]
		if name[0]=='g':
			comment = file.readline().strip().split()[1]
		dimension = file.readline().strip().split()[1]
		edge_weight_type = file.readline().strip().split()[1]
		edge_weight_format = file.readline().strip().split()[1]
		if name[0] != 'g':
			display_data_type = file.readline().strip().split()[1]
		file.readline()
		values = file.readline()
		while values!='EOF\n':
			values = values.split()
			if type(distances)==list:
				distances.append(list(map(int,values)))
			else:
				distances = list(map(int,values))
			
			#if (edge_weight_format=="LOWER_DIAG_ROW"):
			distances = list(flatten(distances))
			values = file.readline()

		while 0 in distances:
				distances.remove(0)
		values = np.array(distances)
		n = int(np.sqrt(len(values)*2))+1
		if edge_weight_format == "LOWER_DIAG_ROW":
			idx = np.tril_indices(n, k=-1, m=n)
		else:
			idx = np.triu_indices(n, k=1, m=n)
			
		matrix = np.zeros((n,n)).astype(int)
		matrix[idx] = values
		neighbours_dict = generate_neighbours(matrix)
		print(matrix)
		print("\n\n")
		first_solution = construct_initial_solution(neighbours_dict)
		print("First solution nodes: ", first_solution)
		cost = calculate_cost(neighbours_dict,first_solution)
		print("Cost solution: ", cost,"\n\n")
		tabu_search(neighbours_dict,first_solution,cost,maxiter,maxcandidate,maxTabu)
		#first_solution = first_solution[:len(first_solution)-(len(first_solution)-maxcandidate)]

def flatten(lis):
    for item in lis:
       if isinstance(item, collections.Iterable) and not isinstance(item, str):
           for x in flatten(item):
               yield x
       else:        
           yield item

def generate_neighbours(matrix):
	dict_of_neighbours = {}
	for i in range(len(matrix)):
		for j in range(len(matrix[i])):
			if(matrix[i][j]!=0):
				if j not in dict_of_neighbours:
					_list = list()
					_list.append([i,matrix[i][j]])
					dict_of_neighbours[j]=_list
				else:
					dict_of_neighbours[j].append([i,matrix[i][j]])
				if i not in dict_of_neighbours:
					_list = list()
					_list.append([j,matrix[i][j]])
					dict_of_neighbours[i]=_list
				else:
					dict_of_neighbours[i].append([i,matrix[i][j]])
	return dict_of_neighbours

def construct_initial_solution(neighbours_dict):
	nodes_visted = list()
	for key in neighbours_dict:
		nodes_visted.append(key)
	size=len(nodes_visted)
	for index in range(size):
		shuffleIndex = random.randrange(index,size)
		nodes_visted[shuffleIndex], nodes_visted[index]= nodes_visted[index], nodes_visted[shuffleIndex]
	return nodes_visted

def calculate_cost(matrix_cost,node_list):
	distance_cost = 0
	for index_node in range(len(node_list)):
		if (index_node < len(node_list)-1):
			if (node_list[index_node]>node_list[index_node+1]):
				index = matrix_cost.get(node_list[index_node+1])
				value = index[node_list[index_node]-1]
				distance_cost+=value[1]
			else:
				index = matrix_cost.get(node_list[index_node])
				value = index[node_list[index_node+1]-1]
				distance_cost+=value[1]
	return distance_cost

def isTabu(perm, tabuList):
    result = False
    size = len(perm)
    for index, edge in enumerate(perm):
        if index == size-1:
            edge2 = perm[0]
        else:
            edge2 = perm[index+1]
        if [edge, edge2] in tabuList:
            result = True
            break
    return result  

def stochasticTwoOptWithEdges(perm):
    result = perm[:] # make a copy
    #print("gere")
    size = len(result)
    # select indices of two random points in the tour
    p1, p2 = random.randrange(0,size), random.randrange(0,size)
    # do this so as not to overshoot tour boundaries
    exclude = set([p1])
    if p1 == 0:
        exclude.add(size-1)
    else:
        exclude.add(p1-1)
    
    if p1 == size-1:
        exclude.add(0)
    else:
        exclude.add(p1+1) 
                       
    while p2 in exclude:
        p2 = random.randrange(0,size)

    # to ensure we always have p1<p2        
    if p2<p1:
        p1, p2 = p2, p1
     
    # now reverse the tour segment between p1 and p2   
    result[p1:p2] = reversed(result[p1:p2])
    
    return result, [[perm[p1-1],perm[p1]],[perm[p2-1],perm[p2]]]


def generate_candidate(best,tabu_list,solution,neighbours_dict):
	permutation, edges, result = None, None, {}
	while permutation == None or isTabu(best["Permutation"], tabu_list):
		permutation, edges = stochasticTwoOptWithEdges(best["Permutation"])
	candidate ={}    
	candidate["Permutation"] = permutation
	test = candidate.get("Permutation")
	candidate["Cost"] = calculate_cost(neighbours_dict,test)
	result["Candidate"] = candidate
	result["Edges"] = edges
	return result

def locateBestCandidate(candidates):
    candidates.sort(key = lambda c: c["Candidate"]["Cost"])
    best, edges = candidates[0]["Candidate"], candidates[0]["Edges"]
    return best, edges 

def tabu_search(neighbours_dict,first_solution,first_cost,maxiter,maxcandidate,maxTabu):
	tabu_list=list()
	best ={}
	best["Permutation"] = first_solution
	best["Cost"] = first_cost
	while maxiter >0:
		print("Iteraction ",maxiter)
		candidates = []
		for index_node in range(0,maxcandidate):
			candidates.append(generate_candidate(best, tabu_list, first_solution,neighbours_dict))
		bestCandidate, bestCandidateEdges = locateBestCandidate(candidates)
		print("\tCandidate: ",bestCandidate.get("Permutation")," Cost: ",bestCandidate.get("Cost"))
		print("\tCurrent Best: ", best.get("Permutation")," Cost: ",best.get("Cost"))
		if bestCandidate["Cost"] < best["Cost"]:
		    # set current to the best, so thatwe can continue iteration
		    print("\tSwitching neighborhood")
		    best = bestCandidate
		    # update tabu list
		    for edge in bestCandidateEdges:
		        if len(tabu_list) < maxTabu:
		            tabu_list.append(edge)
		
		print("\n\n")
		maxiter-=1

	print("Best arrange: ",best)

arguments()