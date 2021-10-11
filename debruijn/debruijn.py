#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
import matplotlib.pyplot as plt
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics

__author__ = "Soula Sarah"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "spequez@gmail.com"
__status__ = "Developpement"

def isfile(path):
	"""Check if path is an existing file.
	  :Parameters:
		  path: Path to the file
	"""
	if not os.path.isfile(path):
		if os.path.isdir(path):
			msg = "{0} is a directory".format(path)
		else:
			msg = "{0} does not exist.".format(path)
		raise argparse.ArgumentTypeError(msg)
	return path


def get_arguments():
	"""Retrieves the arguments of the program.
	  Returns: An object that contains the arguments
	"""
	# Parsing arguments
	parser = argparse.ArgumentParser(description=__doc__, usage=
									 "{0} -h"
									 .format(sys.argv[0]))
	parser.add_argument('-i', dest='fastq_file', type=isfile,
						required=True, help="Fastq file")
	parser.add_argument('-k', dest='kmer_size', type=int,
						default=22, help="K-mer size (default 21)")
	parser.add_argument('-o', dest='output_file', type=str,
						default=os.curdir + os.sep + "contigs.fasta",
						help="Output contigs in fasta file")
	parser.add_argument('-f', dest='graphimg_file', type=str,
						help="Save graph as image (png)")
	return parser.parse_args()

def read_fastq(fastq_file):
	fich = open(fastq_file,"r")
	lignes = fich.readlines()
	for i in range(1,len(lignes),4): 
		yield lignes[i][:-1]
		
	fich.close()

def cut_kmer(read, kmer_size):
	a = 0 
	for i in range(len(read)-kmer_size):
		yield read[i:(i+kmer_size)]
		a+=1
	yield read[a:len(read)]



def build_kmer_dict(fastq_file, kmer_size):
	gene = read_fastq(fastq_file)
	dict= {}
	for j in gene: 
			gene_cut = cut_kmer(j,kmer_size)
			for k in gene_cut:
				if (k in dict):
					dict[k]+=  1    
				else:
					dict[k]=  1
	return dict



def build_graph(kmer_dict):
	G = nx.DiGraph()
	for kmer,val in kmer_dict.items():
		G.add_edge(kmer[0:-1],kmer[1:],weight=val)  
	return G


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):

	if (delete_entry_node ==True and delete_sink_node == True ):
		for i in range(len(path_list)):
			graph.remove_nodes_from(path_list[i])
	if (delete_entry_node ==True and delete_sink_node == False ):
		for i in range(len(path_list)):
			graph.remove_nodes_from(path_list[i][:-1])
	if (delete_entry_node ==False and delete_sink_node == True):
		for i in range(len(path_list)):
			graph.remove_nodes_from(path_list[i][1:])
	if (delete_entry_node ==False and delete_sink_node == False):
		for i in range(len(path_list)):
			graph.remove_nodes_from(path_list[i][1:-1])
	return graph
	
def std(data):
	n = len(data)
	
	moy= statistics.mean(data)
	var=0
	for i in range(n):
		var += abs(data[i]-moy)**2
	var = var/(n-1)
	return var**(0.5)


def select_best_path(graph, path_list, path_length, weight_avg_list, 
					 delete_entry_node=False, delete_sink_node=False):
	if len(weight_avg_list)>=2 and std(weight_avg_list)>0:
		ind = weight_avg_list.index(max(weight_avg_list))
	elif len(path_length)>=2 and std(path_length)>0: 
			ind = path_length.index(max(path_length))
	else: 
			ind = random.randint(0,(len(path_list)-1))
	chemin = path_list[ind]
	path_list.pop(ind)
	graph=remove_paths(graph,path_list,delete_entry_node,delete_sink_node)
	return graph

def path_average_weight(graph, path):
	gene = graph.subgraph(path).edges(data=True)
	return statistics.mean([d['weight'] for (u, v, d) in gene])

def solve_bubble(graph, ancestor_node, descendant_node):
	path_list = list(nx.all_simple_paths(graph,ancestor_node,descendant_node))
	path_length = []
	path_average = []
	for i in range(len(path_list)):
		path_length.append(len(path_list[i]))
		path_average.append(path_average_weight(graph,path_list[i])) 
		graph=select_best_path (graph,path_list,path_length,path_average, False, False)
		return graph


def simplify_bubbles(graph):
	bubble = False 
	gene_nodes = graph.nodes()
	for node in gene_nodes:
		if node in graph: 
			prede = list(graph.predecessors(node))
			if len(prede)>1: 
				for i in range(len(prede)):
					for j in range(len(prede)):
						noeud_ances = nx.lowest_common_ancestor(graph,prede[i] , prede[j])
						if noeud_ances!= None: 
							bubble = True 
							break 
	if bubble:
		graph = simplify_bubbles(solve_bubble(graph,noeud_ances, node))

def solve_entry_tips(graph, starting_nodes):
	liste_noeud = graph.nodes()
	for noeud in liste_noeud: 
		if noeud in graph:
			prede = graph.predecessors(noeud)
			liste_pred_entr=[]
			for pred in prede : 
				if pred in starting_nodes:
					liste_pred_entr.append(pred)
			if len(liste_pred_entr)== 2: 
				liste_pred_entr.insert(1,noeud)
				graph.remove_nodes_from(liste_pred_entr)
				starting_nodes.remove(liste_pred_entr[0])
				starting_nodes.remove(liste_pred_entr[2])
				solve_entry_tips(graph,starting_nodes)
			
def solve_out_tips(graph, ending_nodes):
	liste_noeud = graph.nodes()
	for noeud in liste_noeud: 
		if noeud in graph:
			prede = graph.successors(noeud)
			liste_pred_entr=[]
			for pred in prede : 
				if pred in ending_nodes:
					liste_pred_entr.append(pred)
			if len(liste_pred_entr)== 2: 
				liste_pred_entr.insert(1,noeud)
				graph.remove_nodes_from(liste_pred_entr)
				ending_nodes.remove(liste_pred_entr[0])
				ending_nodes.remove(liste_pred_entr[2])
				solve_out_tips(graph,ending_nodes)

def get_starting_nodes(graph):
	liste_noeud = graph.nodes()
	liste_start = []
	for el in liste_noeud: 
		gene = graph.predecessors(el)
		if not any(gene): 
			liste_start.append(el)
	return liste_start


def get_sink_nodes(graph):
	liste_noeud = graph.nodes()
	liste_end = []
	for el in liste_noeud: 
		gene = graph.successors(el)
		if not any(gene): 
			liste_end.append(el)
	return liste_end

def get_contigs(graph, starting_nodes, ending_nodes):
	liste_contig = []
	kmer=len(starting_nodes[0])
	for i in range(len(starting_nodes)):
		for j in range(len(ending_nodes)):
			if nx.has_path(graph,starting_nodes[i],ending_nodes[j]):
				gene = nx.all_simple_paths(graph,starting_nodes[i],ending_nodes[j])
				for el in gene:
					chaine=""
					for k in range(len(el)):
						if k%kmer==0:
							chaine+=el[k]				
					liste_contig.append(( chaine,len(chaine)))
	return liste_contig
				

def save_contigs(contigs_list, output_file):
	with open(output_file,"w") as filout:
		for i in range(len(contigs_list)): 
			filout.write(f">contig_{i} len={contigs_list[i][1]}\n")
			filout.write(fill(contigs_list[i][0], width=80))
			filout.write("\n")



def fill(text, width=80):
	"""Split text with a line return to respect fasta format"""
	return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def draw_graph(graph, graphimg_file):
	"""Draw the graph
	"""                                    
	fig, ax = plt.subplots()
	elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
	#print(elarge)
	esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
	#print(elarge)
	# Draw the graph with networkx
	#pos=nx.spring_layout(graph)
	pos = nx.random_layout(graph)
	nx.draw_networkx_nodes(graph, pos, node_size=6)
	nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
	nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
						   edge_color='b', style='dashed')
	#nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
	# save image
	plt.savefig(graphimg_file)


def save_graph(graph, graph_file):
	"""Save the graph with pickle
	"""
	with open(graph_file, "wt") as save:
			pickle.dump(graph, save)


#==============================================================
# Main program
#==============================================================
def main():
	"""
	Main program function
	"""
	# Get arguments
	#args = get_arguments()
	dic = build_kmer_dict("data/eva71_two_reads.fq",6)
	graph = build_graph(dic)
	#draw_graph(graph,"graph.png")
	start_node = get_starting_nodes(graph)
	end_node = get_sink_nodes(graph)
	#print(start_node)
	#print(end_node)
	list_contig = get_contigs(graph, start_node,end_node)
	#save_contigs(list_contig,"contigs.txt")
	#print()
   

	# Fonctions de dessin du graphe
	# A decommenter si vous souhaitez visualiser un petit 
	# graphe
	# Plot the graph
	# if args.graphimg_file:
	#     draw_graph(graph, args.graphimg_file)
	# Save the graph in file
	# if args.graph_file:
	#     save_graph(graph, args.graph_file)


if __name__ == '__main__':
	main()
