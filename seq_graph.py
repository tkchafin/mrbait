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
    return(nx.maximal_independent_set(G, nodes=None))


#Function
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


def plotColorNodes(G, listnodes):
    color_map = []
    for node in G:
        if node in listnodes:
            color_map.append("red")
        else:
            color_map.append("black")
    nx.draw(G, node_color = color_map, with_labels=True)
    plt.show()

#Tests of functions
example_10 = [(1,2),(2,4),(1,3),(1,7),(3,2),(1,4),(5,6),(6,8),(3,7),(4,8),(9,10)]
example_100 = [(19,29),(28,48),(17,36),(16,72),(33,2),(1,47),(55,66),(62,87),(53,57),(64,68),(9,100),
(11,22),(24,46),(11,32),(89,78),(31,24),(19,45),(54,6),(16,88),(3,7),(4,88),(95,43),
(11,28),(27,4),(1,38),(13,7),(3,2),(1,48),(49,57),(61,8),(98,79),(81,80),(97,100),
(12,29),(26,4),(1,37),(1,71),(39,2),(1,47),(50,58),(36,8),(63,78),(24,82),(96,100),
(13,30),(25,4),(78,36),(12,7),(40,2),(1,46),(56,59),(61,99),(3,77),(4,83),(95,11),
(14,12),(24,4),(1,35),(14,15),(3,2),(1,42),(55,60),(6,100),(3,76),(4,84),(92,94),
(15,2),(23,4),(2,31),(1,71),(3,2),(1,43),(51,6),(63,64),(70,7),(4,85),(90,93),
(16,23),(21,34),(14,32),(12,7),(12,13),(1,41),(52,61),(62,8),(71,72),(4,86),(91,10),
(17,21),(22,64),(27,33),(14,7),(83,72),(1,45),(53,69),(65,8),(74,73),(4,87),(89,10),
(18,22),(20,4),(59,34),(1,45),(91,75),(19,44),(54,67),(66,68),(31,75),(45,18),(90,10)
]

G10 = multiGraphFromList(example_10)
G100 = multiGraphFromList(example_100)

z = nx.utils.create_degree_sequence(20,powerlaw_sequence)
G = nx.configuration_model(z)
G=nx.Graph(G)
G.remove_edges_from(G.selfloop_edges())

#approximateIndependentSet(100,100)
C = naiveIndependentSet(G)

nodes = C.nodes()
#nodes = approximateIndependentSet(G)
#plotColorNodes(G, nodes)
