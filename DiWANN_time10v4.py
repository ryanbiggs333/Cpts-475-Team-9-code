#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
from igraph import Graph, mean
from igraph import VertexSeq
from igraph import EdgeSeq
from igraph import summary
from igraph import plot
from igraph import GraphBase
import igraph as ig
from igraph import VertexClustering
from igraph import clustering
import time
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap
import matplotlib.cm as cm
from sklearn.cluster import AgglomerativeClustering
import scipy.cluster.hierarchy as sch
#import louvain

#This function returns a sanitized protein with spaces and digits removed, and case set to lower		
def sanitize(string):
	string=re.sub(r'[^\w]','',string)	
	string=re.sub(r'[\d_]','',string)		
	string=string.lower()
	return string

#returns a list as a string with sep separators.
def listtostring(list, sep):
	string=""
	if len(list)>0:
		string+=list[0]#.decode('utf-8')
	for item in list[1:]:
		string+=sep+item#.decode('utf-8')
	return string

#function to read fasta files
def readFasta(file, keepFile="NA"):
	import re
	sequences=[]
	kp=set()
	with open(file,'r') as input:
		if keepFile!="NA":
			with open(keepFile,'r') as keep:
				for line in keep:
					kp.add(line[:-1])
		#print kp
		for line in input:
			if line[0]!=">":
				if len(kp)>0:
					if sequences[-1][0] in kp:
						#print "keep"
						sequences[-1][1]+=sanitize(line)
					else:
						#print sequences[-1], "not on keep list"
						continue
				else:
					sequences[-1][1]+=sanitize(line)
			else:
				if len(sequences)>0 and sequences[-1][1]=="":
					del sequences[-1]
				sequences.append( [re.sub('[^0-9]', '',line[1:]),""] )
	print(len(sequences), sequences[0])
	return sequences

#Edit Distance
def levenshteinDistance(s1,s2):
		
    if len(s1) > len(s2):
        s1,s2 = s2,s1
    distances = range(len(s1) + 1)
    for index2,char2 in enumerate(s2):
        newDistances = [index2+1]
        for index1,char1 in enumerate(s1):
            #scoremod=matrix[(char1,char2)]-4 #Value between 0 and -8
            if char1 == char2:
                newDistances.append(distances[index1])
                #newDistances.append(distances[index1]-scoremod)
            else:
                newDistances.append(1 + min((distances[index1],
                                             distances[index1+1],
                                             newDistances[-1])))
                #newDistances.append(1 + min((distances[index1],
                #                             distances[index1+1],
                #                             newDistances[-1])) - scoremod)
        distances = newDistances
#   print distances[-1]
    return distances[-1]
#END

#bounded edit distance    
def levenshteinDistanceThresh(s1,s2,k):
#    print "LD:", k
#    print s1
#    print s2

    from Bio.SubsMat import MatrixInfo as matlist
    matrix = matlist.blosum62

    if len(s1) > len(s2):
        s1,s2 = s2,s1
    distances = range(len(s1) + 1)
    for index2,char2 in enumerate(s2):
        newDistances = [index2+1]
        for index1,char1 in enumerate(s1):
            if abs(index1-index2) > k:
                newDistances.append(sys.maxsize)
            elif char1 == char2:
                newDistances.append(distances[index1])
            else:
                newDistances.append(1 + min((distances[index1],distances[index1+1],newDistances[-1])))	
        distances = newDistances

        if min(distances)>k:
            return sys.maxsize		
    return distances[-1]

#DiWANN graph construction
def makeGraphEfficiently(species, repeatmatrix, filename, sparse=False):
	
	G = Graph(directed=True)
	G.add_vertices(len(species))
	
	#summary (G)
	
	names=[]
	for r in labelss:
		names.append('; '.join(r))
	#G.vs["label"] = names
	
###
	for row in range(len(species)):
		#populate graph
		#print len(repeatmatrix[row])
		if sparse==True:
			import scipy
			m=min(list(scipy.sparse.find(repeatmatrix.tocsc()[:,row])[2]) + list(scipy.sparse.find(repeatmatrix.tocsc()[row,:])[2]))
			edgeTo = getSparseMinED(repeatmatrix, row, m)
			#return
		else:
			m=min(l for l in repeatmatrix[row] if l > 0)
			edgeTo = getMinED(repeatmatrix, row, m)
		#print edgeTo
		for e in edgeTo:
			G.add_edge(e[0], e[1], weight=m)
		
	#summary(G)
	G.simplify(combine_edges=min)
	summary(G)
	G.write(filename+"DWRNgraphEff.gml","gml")
	#summary(G)
	return G

#threshold based network construction
def makeGraph(species, RepeatMatrix, filename):

	G_new = Graph(directed=False)
	G_new.add_vertices(len(species))
	
	G = Graph()
	G.add_vertices(len(species))
	
	names=[]
	for r in labelss:
		names.append('; '.join(r))

	#G_new.vs["label"] = names
	#G.vs["label"] = names
	
	its=len(species)*len(species)
	it=0
	for i in range(len(species)):
		# print RepeatMatrix[i]
		minED = min(l for l in RepeatMatrix[i] if l > 0)
		
		for j in range(len(species)):
			#print "adding"
			if (RepeatMatrix[i][j] == minED):
				G_new.add_edge(i,j,weight=RepeatMatrix[i][j])
				#weights.append(RepeatMatrix[i][j])
			if i>j and RepeatMatrix[i][j]<2000000:
				G.add_edge(i,j,weight=RepeatMatrix[i][j])
			#print "added"
			it+=1
		print(round((it*1.0)/its,3))
	summary(G_new)
	summary(G)
	
	
	G_new.write(filename+"DWN-base.gml","gml")
	G.write(filename+"Thresh.gml","gml")
 
	return G

#finds lowest in a vector
def getMinED(repeatmatrix, row, minVal):
#	minVal=sys.maxint
	
	#print repeatmatrix
	ret=[]
	
	for col in range(len(repeatmatrix[row])):
#		#print "a"
#		if repeatmatrix[row][col] < minVal:
#			minVal = repeatmatrix[row][col]
#			ret = [(row,col)]
#		elif repeatmatrix[row][col] == minVal:
		if repeatmatrix[row][col] == minVal:
			ret.append((row,col))
		if repeatmatrix[row][col] < minVal and repeatmatrix[row][col] < 0:
			print("Warning, min passed was wrong")
			ret = [(row,col)]
			minVal=repeatmatrix[row][col] 
	
	#print minVal
	#print ret
	return ret
	
def getSparseMinED(repeatmatrix, row, minVal):
#	minVal=sys.maxint
	
	#print repeatmatrix
	ret=[]
	
	for col in range(len(repeatmatrix[row])):
#		#print "a"
#		if repeatmatrix[row][col] < minVal:
#			minVal = repeatmatrix[row][col]
#			ret = [(row,col)]
#		elif repeatmatrix[row][col] == minVal:
		if repeatmatrix[row,col] == minVal:
			ret.append((row,col))
		if repeatmatrix[row,col] < minVal and repeatmatrix[row,col] < 0:
			print("Warning, min passed was wrong")
			ret = [(row,col)]
			minVal=repeatmatrix[row,col] 
	
	#print minVal
	#print ret
	return ret

#creats DiWANN similarity matrix
def createMinMatrix(species, useBoundedLD=False):
    #print("started creatMinMatrix")
    from Bio.SubsMat import MatrixInfo as matlist
    from Bio import pairwise2
    matrix = matlist.blosum62
	#print(matrix)
    skippedCells=0	
    repeatmatrix=[]
    its=(len(species)*len(species))/2
    it=len(species)
    #print("its=",its)
    #print("it=",it)
    for row in range(len(species)):
        #Add new row to matrix
        repeatmatrix.append([])
        repeatmatrix[row].extend([sys.maxsize]*len(species))
    #print(repeatmatrix)
    for row in range(len(species)):
        #print "row:", row
        if row>0:
            print("Percent complete:", round((it*1.0)/its,3))
            #calculate ranges
            maxED=[]
            minED=[]
            if row>1:
                rowMin=min(repeatmatrix[row][0:row-1])
                #print repeatmatrix[row][0:row-1]
                #print rowMin
            else:
                rowMin=sys.maxsize
			
            for col in range(row+1,len(species)):
                minED.append(abs(repeatmatrix[0][col]-repeatmatrix[0][row]))
                maxED.append(repeatmatrix[0][col]+repeatmatrix[0][row])
			#then only compare possible mins
			#print row, len(minED)
			
            if len(maxED)>0:
                lowestMax = min(maxED)
            else:
                lowestMax = sys.maxsize
                #print lowestMax
			
            for col in range(row+1,len(species)):
                it+=1
                colMin = min(repeatmatrix[col])
                #print colMin
                #if col - (row+1) < 0:
                #print "ERROR"
                if (minED[col-(row+1)] > lowestMax or minED[col-(row+1)] > rowMin) and (minED[col-(row+1)] > colMin): 
                    repeatmatrix[row][col]=sys.maxsize
                    repeatmatrix[col][row]=sys.maxsize
					#print "skipping ", row, col
                    skippedCells+=1
                else:
                    if useBoundedLD == True:
                        repeatmatrix[row][col]=levenshteinDistanceThresh(species[row][0],species[col][0], max(colMin,rowMin))
					 #repeatmatrix[row][col]=pairwise2.align.globaldx(species.repeats[row].sequence.upper(), species.repeats[col].sequence.upper(), matrix, score_only=True, one_alignment_only=True)
                        if repeatmatrix[row][col] > max(colMin,rowMin):
                            repeatmatrix[row][col]=sys.maxsize
                        #print max(colMin,rowMin)
                    else:
                        repeatmatrix[row][col]=levenshteinDistance(species[row][0],species[col][0])
                        #repeatmatrix[row][col]=pairwise2.align.globaldx(species.repeats[row].sequence.upper(), species.repeats[col].sequence.upper(), matrix, score_only=True, one_alignment_only=True)
                    repeatmatrix[col][row]=repeatmatrix[row][col]
                    #if repeatmatrix[row][col] == 0:
                        #print "WARNING!!!"
                    if repeatmatrix[row][col] < rowMin:
                        rowMin=repeatmatrix[row][col]

        else:
            for col in range(len(species)):
                repeatmatrix[row][col]=levenshteinDistance(species[row][0],species[col][0])
                #repeatmatrix[row][col]=pairwise2.align.globaldx(species.repeats[row].sequence.upper(), species.repeats[col].sequence.upper(), matrix, score_only=True, one_alignment_only=True)
                #print repeatmatrix[col][row]
                repeatmatrix[col][row]=repeatmatrix[row][col]
            repeatmatrix[0][0] = sys.maxsize

    # print("Skipped computing", skippedCells, "of", ((len(species)*len(species))-len(species))/2, "cells")
    return repeatmatrix

#Creats similarity matric for threshold-based methods
def createRepeatMatrix(species, thresh=2):
    from Bio.SubsMat import MatrixInfo as matlist
    from Bio import pairwise2
    matrix = matlist.blosum62

    repeatmatrix=[]
    for row in range(len(species)):
        #Add new row to matrix
        repeatmatrix.append([ listtostring(labelss[row],";") ])
        repeatmatrix[row].extend([0]*len(species))
    its=len(species)*len(species)/2
    i=0
    for row in range(len(species)):
        for col in range(row,len(species)):
            #repeatmatrix[row][col]=levenshteinDistance(species[row][0],species[col][0])
            repeatmatrix[row][col]=levenshteinDistanceThresh(species[row][0],species[col][0],thresh)
            #repeatmatrix[row][col]=pairwise2.align.globaldx(species[row][0].upper(), species[col][0].upper(), matrix, score_only=True, one_alignment_only=True)
            repeatmatrix[col][row]=repeatmatrix[row][col]
            
            i+=1
        print("Percent complete:", round((i*1.0)/its,3))

    return repeatmatrix

def makeNetwork(filename,thresh,bthresh,type="", sequences=None,labels = None):
    graph=""
    if type=="PB":
        print("Prune + bound*****************")
        print(filename)
        start = time.time()
        C = createMinMatrix(sequences, False)
        C = np.asarray(C)
        end = time.time()
       # y_hc=create_dendrogram(C)
        print("Pruning and Bounding:", (end-start))
        #MatrixtoCSV(C,"DWRN"+filename+".csv")
        graph = makeGraphEfficiently(sequences, C, filename+"minNet")
        summary(graph)    
    return graph

style = {}
style["edge_curved"] = True
style["margin"] = 50
style["edge_arrow_size"]=0.1
style["vertex.label.cex"]=0.4
style["vertex.label.family"]="Helvetica"
style["vetrex.label.color"]='black'
style["edge_width"]=0.5
style["edge_color"]="black"
style["vertex_size"]=9
style["vertex.label.font"]=1.5


# load csv file
donor_sequence_df = pd.read_csv('donor_art_seq_468v5.csv')

# Check for missing values and drop them
donor_sequence_df = donor_sequence_df.dropna(subset=['seq'])

# Convert 'seq' column to string
donor_sequence_df['seq'] = donor_sequence_df['seq'].astype(str)

# extract sequences and labels
sequences = donor_sequence_df['seq'].apply(sanitize).tolist()
donor_ids = donor_sequence_df['sample id'].tolist()

rfile2 = "network_output"
labelss = donor_ids

#create DiWANN network - our base SNN
diwann_snn = Graph()
diwann_snn = makeNetwork(rfile2,thresh="",bthresh="", type="PB", sequences=sequences, labels=donor_ids)
plot(diwann_snn, 'diwann_network_output.png', **style)

#----------------- tissue visualization -----------------
tissue_types = donor_sequence_df['Project_Type'].tolist()

# Add tissue type as an attribute to each node in the SSN
for i in range(len(diwann_snn.vs)):
    diwann_snn.vs[i]["Project_Type"] = tissue_types[i]

# Create a color map for tissue types
unique_tissues = list(set(tissue_types))
colors = plt.cm.get_cmap('viridis', len(unique_tissues))
tissue_color_map = {tissue: colors(i) for i, tissue in enumerate(unique_tissues)}

# Apply colors to nodes based on tissue type
node_colors = [tissue_color_map[tissue] for tissue in tissue_types]
style["vertex_color"] = node_colors

# output SNN
plot(diwann_snn, 'diwann_network_tissue_overlay.png', **style)

fig, ax = plt.subplots()

# Create a list of handles for the legend, one for each tissue type
handles = [plt.Line2D([0], [0], marker='o', color='w', label=tissue,
                      markerfacecolor=tissue_color_map[tissue], markersize=10)
           for tissue in unique_tissues]

# Create the legend
legend = ax.legend(handles=handles, loc='center', title='Tissue Types')

ax.axis('off')
fig.canvas.draw()
bbox = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())

# Save the legend to a file
plt.savefig('tissue_legend.png', bbox_inches=bbox, pad_inches=0)

#----------------- cluster visualization -----------------
undirected_diwann_snn = diwann_snn.as_undirected()

# Community detection to find clusters on undirected graph
clusters = undirected_diwann_snn.community_multilevel()

diwann_snn.vs['cluster'] = clusters.membership
modularity_score = undirected_diwann_snn.modularity(clusters.membership)
print("Modularity Score:", modularity_score)

# Analyze the composition of the clusters
for cluster_id in set(clusters.membership):
    print(f"Cluster {cluster_id}:")
    cluster_nodes = diwann_snn.vs.select(lambda vertex: vertex['cluster'] == cluster_id)
    tissues_in_cluster = [v['Project_Type'] for v in cluster_nodes]
    print(f"Tissue types in cluster: {set(tissues_in_cluster)}")

# Create a color map for cluster visualization
n_clusters = len(set(clusters.membership))
colors = plt.get_cmap('Accent', n_clusters)

# Print out the color
for cluster_id in set(clusters.membership):
    color = colors(cluster_id)
    print(f"Cluster {cluster_id} color: {color}")

# Define style settings for the graph plotting
style = {
    "edge_curved": False,
    "margin": 50,
    "edge_arrow_size": 0.1,
    "vertex_label_cex": 0.4,
    "vertex_label_family": "Helvetica",
    "vertex_label_color": 'black',
    "edge_width": 0.5,
    "edge_color": "black",
    "vertex_size": 9,
    "vertex_label_font": 1.5,
    "vertex_color": [colors(c) for c in clusters.membership]
}

ig.plot(diwann_snn, 'diwann_network_clusters.png', **style)

fig2, ax2 = plt.subplots()
# Create a list of handles for the legend
handles = [plt.Line2D([0], [0], marker='o', color='w', label=f'Cluster {i}',
                      markerfacecolor=colors(i), markersize=10) for i in range(n_clusters)]

# Create the legend
legend = ax2.legend(handles=handles, loc='upper left', title='Clusters')

ax2.axis('off')
fig2.savefig('legend.png', bbox_inches=legend.get_window_extent().transformed(fig2.dpi_scale_trans.inverted()), pad_inches=0)

#----------------- tumor visualization -----------------

donor_sequence_df['tumor_gene_uni'] = donor_sequence_df['tumor_gene_uni'].apply(lambda x: x.split(','))

# list of unique genes
all_genes2 = [gene for sublist in donor_sequence_df['tumor_gene_uni'] for gene in sublist]
unique_genes2 = set(all_genes2)

# Create a color map for genes
gene_to_color2 = {gene: plt.cm.tab20(i % 20) for i, gene in enumerate(unique_genes2)}
node_colors_by_gene2 = [gene_to_color2[gene_list[0]] if gene_list else (0, 0, 0, 1)  # Default to black if no gene
                       for gene_list in donor_sequence_df['tumor_gene_uni']]

style = {
    "vertex_color": node_colors_by_gene2
}

# Plot the network with gene-based color coding
ig.plot(diwann_snn, 'diwann_network_with_genes2.png', **style)

def save_gene_legend(gene_to_color, filename):
    handles = [Line2D([0], [0], marker='o', color=gene_to_color[gene], label=gene, markersize=10)
               for gene in gene_to_color]
    fig, ax = plt.subplots()
    legend = ax.legend(handles=handles, loc='center', title='Gene Tumors')
    ax.axis('off')
    fig.canvas.draw()
    bbox = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    
    # Save the legend to a file with the bounding box as the dimensions
    fig.savefig(filename, bbox_inches=bbox, pad_inches=0)
    plt.close(fig)

save_gene_legend(gene_to_color2, 'gene_legend2.png')

#----------------- cancer type visualization -----------------
cancer_types = donor_sequence_df['Cancer_Type'].tolist()

# Add cancer type as an attribute to each node in the SSN
for i in range(len(diwann_snn.vs)):
    diwann_snn.vs[i]["Cancer_Type"] = cancer_types[i]

# Create a color map for cancer types
unique_cancers = list(set(cancer_types))
colors = plt.cm.get_cmap('viridis', len(unique_cancers))
cancer_color_map = {cancer: colors(i) for i, cancer in enumerate(unique_cancers)}

# Apply colors to nodes based on cancer type
node_colors = [cancer_color_map[cancer] for cancer in cancer_types]
style["vertex_color"] = node_colors

# output SNN
plot(diwann_snn, 'diwann_network_cancer_overlay.png', **style)

fig, ax = plt.subplots()

# Create a list of handles for the legend, one for each cancer type
handles = [plt.Line2D([0], [0], marker='o', color='w', label=cancer,
                      markerfacecolor=cancer_color_map[cancer], markersize=10)
           for cancer in unique_cancers]

# Create the legend
legend = ax.legend(handles=handles, loc='center', title='Cancer Types')

ax.axis('off')
fig.canvas.draw()
bbox = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())

# Save the legend to a file
plt.savefig('cancer_legend.png', bbox_inches=bbox, pad_inches=0)


#Visualization example using a list color_mapT10 containing node colors according to the state label
'''
bstyle = {}
bstyle["edge_curved"] = False
bstyle["margin"] = 50
bstyle["edge_arrow_size"]=0.1
bstyle["vertex.label.cex"]=0.2
bstyle["vertex.label.family"]="Helvetica"
bstyle["vetrex.label.color"]='black'
bstyle["edge_width"]=0.2
#bstyle["edge_color"]=black
size=[]
#seq_label=[]
for i in range(len(uni_count_T10)):
    if uni_count_T10[i] == 1:
        size.append(7)
    if uni_count_T10[i] > 1 and uni_count_T10[i] <= 10:
        size.append(7)
    if uni_count_T10[i] > 10 and uni_count_T10[i] <= 60:
        size.append(8)
    if uni_count_T10[i] > 60 and uni_count_T10[i] <= 500:
        size.append(10)
    if uni_count_T10[i] > 500:
        size.append(14)
        #seq_label.append(uni_labels_AllT[i])
shape=[]
for a in aa10_95:
    if a == 'I':
        shape.append('triangle-down')
    if a == 'T':
        shape.append('circle')
    if a == 'X':
        shape.append('rectangle')
bstyle["vertex_shape"] = shape
bstyle["vertex_size"]=size
bstyle["vertex.label.font"]=0.9
bstyle["vertex_color"] = color_mapT10
'''