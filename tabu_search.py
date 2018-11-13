import argparse, sys, platform, os, operator, collections
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
	parser.add_argument('--maxlist','-ml',type=int,required=True,help='Tamanho maximo da lista')
	parser.add_argument('--maxcandidate','-mc',type=int,required=True,help='Quantidade maxima de candidados')
	args = parser.parse_args()
	global path
	if platform.system() == "Windows":
		path = path +r'\\'+args.file
	elif platform.system() == "Linux":
		path = path  + '/' + args.file
	decode_file()

def decode_file():
	global path, name, comment, dimension, edge_weight_type, edge_weight_format,distances
	with open(path,'r') as file:
		name = file.readline().strip().split()[1]
		type_file = file.readline().strip().split()[1]
		comment = file.readline().strip().split()[1]
		dimension = file.readline().strip().split()[1]
		edge_weight_type = file.readline().strip().split()[1]
		edge_weight_format = file.readline().strip().split()[1]
		if name[0] != 'g':
			display_data_type = file.readline().strip().split()[1]
			print(display_data_type)
		file.readline()
		values = file.readline()
		while values!='EOF\n':
			values = values.split()
			if type(distances)==list:
				distances.append(list(map(int,values)))
			else:
				distances = list(map(int,values))
			
			if (edge_weight_format=="LOWER_DIAG_ROW"):
				distances = list(flatten(distances))
			values = file.readline()

		if edge_weight_format == "LOWER_DIAG_ROW":
			while 0 in distances:
				distances.remove(0)
			values = np.array(distances)
			n = int(np.sqrt(len(values)*2))+1
			idx = np.tril_indices(n, k=-1, m=n)
			matrix = np.zeros((n,n)).astype(int)
			matrix[idx] = values
			print(matrix)
			print("\n\n")
			neighbours_dict = generate_neighbours(matrix)

		first_solution,distance_of_first_solution = generate_first_solution(neighbours_dict,1)
		print("First solution\ndistance:",distance_of_first_solution,"\n")
		print(first_solution)

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

def generate_first_solution(neighbours_dict,node):
	start_node = node
	end_node = start_node
	first_solution = []
	visiting = start_node
	distance_of_first_solution = 0

	while visiting not in first_solution:
		minim = 10000
		for k in neighbours_dict[visiting]:
			if k[1] < minim and k[0] not in first_solution:
				minim = k[1]
				best_node=k[0]
		first_solution.append(visiting)
		distance_of_first_solution+=minim
		visiting = best_node
	first_solution.append(end_node)
	position=0
	for k in neighbours_dict[first_solution[-1]]:
		if k[0]==start_node:
			break

	distance_of_first_solution += abs(neighbours_dict[first_solution[-2]][position][1]-10000)
	return first_solution,distance_of_first_solution

def find_neighborhood(solution,neighbours_dict):
	neighbours_of_solution = []
	for n in solution[1:-1]:
		idx1 = solution.index(n)
		for kn in solution[1:-1]:
			idx2 = solution.index(kn)
			if n == kn:
				continue
			_tmp = copy.deepcopy(solution)
			_tmp[idx1]=kn
			_tmp[idx2]=n
			distance=0
			for k in _tmp[:-1]:
				next_node = _tmp[_tmp.index(k)+1]
				for i in neighbours_dict[k]:
					if i[0]==next_node:
						distance+=i[1]
			_tmp.append(distance)
			if _tmp not in neighbours_of_solution:
				neighbours_of_solution.append(_tmp)
	index_of_last_item_in_the_list=len(neighbours_of_solution[0])-1
	neighbours_of_solution.sort(key=lambda x:x [index_of_last_item_in_the_list])

arguments()