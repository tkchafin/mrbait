#!/usr/bin/python

import networkx as nx
import networkx.algorithms.approximation as nxaa
import matplotlib.pyplot as plt
import numpy as np
from networkx.utils import powerlaw_sequence

"""Code for ATTEMPTING to approximate the maximal independent set in a graph
of conflicting sequences (e.g. aligned > threshold in pairwise alignment).
Unfortunately, this problem is NP-hard and can't be done efficiently... """

"""Conclusions: My naive version seems to be faster somehow."""

#TODO: Test that nodes retained in NAIVE definitely dont have any edges!!!

def multiGraphFromList(data):
	G = nx.MultiGraph()
	G.add_edges_from(data)
	return(G)

#Function to use the built-in independent set function in networkx
def approximateIndependentSet(G):
	kept = nx.maximal_independent_set(G, nodes=None)
	return(kept)

#Function that loops through all edges, and resolves by keeping the node with the least neighbors
#If there is a conflict, we keep the right node of each edge by default
def naiveIndependentSet(G):
	#Make a copy of graph
	C = G.copy()
	#Loop through ALL edges
	for n in G.edges_iter():
		#If either node has been trimmed from Copy, skip.
		if C.has_edge(n[0],n[1]):
			right = n[1]
			left = n[0]
			right_n = len(C.neighbors(right))
			left_n = len(C.neighbors(left))
			#print("Right neighbor <",right,"> has ", right_n, " connections.")
			#print("Left neighbor <",left,"> has ", left_n, " connections.")
			#Remove right if it has more neighbors, otherwise remove left
			if (right_n > left_n):
				C.remove_node(right)
			else:
				C.remove_node(left)
	return(C)

#Function that loops through all edges, and resolves by keeping the node with the least neighbors
#If there is a conflict, we keep the node which has the most SNPs in proximity (rom weights df)
def weightedNaiveIndependentSet(G, weights):
	#Make a copy of graph
	C = G.copy()
	print(weights)
	#Loop through ALL edges
	for n in G.edges_iter():
		#If either node has been trimmed from Copy, skip.
		if C.has_edge(n[0],n[1]):
			right = n[1]
			left = n[0]
			right_n = len(C.neighbors(right))
			left_n = len(C.neighbors(left))
			#print("Right neighbor <",right,"> has ", right_n, " connections.")
			#print("Left neighbor <",left,"> has ", left_n, " connections.")
			#Remove right if it has more neighbors, otherwise remove left
			if (right_n > left_n): #Remove right if left has more neightbors
				C.remove_node(right)
			elif (left_n > right_n): #remove left if right has more neightbors
				C.remove_node(left)
			else: #else they are equal, attempt weighted resolve
				# print("Conflict, getting weights from: ",right, " and ", left)
				# print("Right weight is: ", weights[int(right)])
				# print("Left weight is: ", weights[int(left)])
				right_weight = weights[int(right)]
				left_weight = weights[int(left)]
				if (left_weight > right_weight): #if left has higher weight, remove right
					C.remove_node(right)
				else:
					C.remove_node(left) #otherwise default to remove left
	return(C.nodes()) #returns list of nodes

#Function to plot a complete graph, coloring a list of 'chosen' or 'excluded' (subset) nodes
def plotColorNodes(G, listnodes):
	color_map = []
	for node in G:
		if node in listnodes:
			color_map.append("red")
		else:
			color_map.append("black")
	nx.draw(G, node_color = color_map, with_labels=True)
	plt.show()

#Function to return a revised blacklist after unweighted edge resolution
#Must be passed a list of edges
def edgeResolveApproximate(edges):
	G = multiGraphFromList(edges) #generate graph from list of edges
	kept = approximateIndependentSet(G) #get retained nodes from independent set construction
	#print("Keeping:", kept)
	remove = []
	blacklist = listFromEdges(edges) #Get list of all "bad" nodes
	#For each bad node, if it wasn't kept: add to REMOVE list so it will be deleted
	for i in blacklist:
		if i not in kept:
			remove.append(i)
			print(i)
	#print("Removing:",remove)
	return(remove) #return list of nodes that were removed.

#Function to return a revised blacklist after weighted edge resolution
#Must be passed a list of edges and a pandas DF of weights
def edgeResolveWeighted(edges, weights):
	G = multiGraphFromList(edges) #generate graph from list of edges
	kept = weightedNaiveIndependentSet(G, weights) #get retained nodes from independent set construction
	#print("Keeping:", kept)
	remove = []
	blacklist = listFromEdges(edges) #Get list of all "bad" nodes
	#For each bad node, if it wasn't kept: add to REMOVE list so it will be deleted
	for i in blacklist:
		if i not in kept:
			remove.append(i)
			print(i)
	#print("Removing:",remove)
	return(remove) #return list of nodes that were removed.


#Function to get a list of member nodes from a list of edges
def listFromEdges(edges):
	ret = [] #empty list
	for i in edges:
		for j in i:
			if j not in ret:
				ret.append(j)
	return(ret)
