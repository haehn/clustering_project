#! /usr/bin/env python

#2008-10-03 exception introduced when too many internal splits, but still error when decomposition is applied

# own libraries
from utils import *
from Graph import *

# external libraries
from string import *
from copy import *
import sys
import math
import time


#================================================================================#  
# Methods for parsing newick trees
#================================================================================#  

def get_splits(newick_string,term=False): #if term: append also terminal splits and branch lengths    
    return splits_for_tree(parse_newick(newick_string),term)

def get_splits_with_decomposition(newick1,newick2,term=False): #if term: append also terminal splits and branch lengths

    
    root1=parse_newick(newick1)
    root2=parse_newick(newick2)

    
    taxad,decomp=root1.decomposition(root2) #returns taxa-dictionary of dummy taxa and decompositions as pairs of root nodes

    all_splits=[] #stores pairs of subtrees to be compared

    for droot1,droot2 in decomp:
        act_splits=[]

        for droot in droot1,droot2:
            act_splits.append(list(splits_for_tree(droot,term)))
        all_splits.append(act_splits)
                
    return taxad,all_splits

def get_split_representation(splits1,splits2): #extract splits only in one tree and compute adjacency matrix
    
    def is_compatible(s1,s2):
        s1r,s1l=s1.split('|') #right and left side
        s2r,s2l=s2.split('|')
        for i in s1r,s1l:
            act1=set(i.split('*'))
            for j in s2r,s2l:
                act2=set(j.split('*')) 
                if len(act1.intersection(act2))==0:return True
        return False

    splits1s=set(splits1).difference(set(splits2)) #splits only in T1
    splits2s=set(splits2).difference(set(splits1)) #splits only in T2

    splits=list(splits1s)+list(splits2s) #now the splits are only those which occur in Exactly one tree, the first dim1 are from T1 and the last dim2 are from T2
    
    dim1=len(splits1s)
    dim2=len(splits2s)
    
    adj=[]
    for i in range(0,dim1):
        adj.append([False]*dim2)
        
    for i in range(0,dim1):
        for j in range(0,dim2):
            c=is_compatible(splits[i],splits[dim1+j])
            adj[i][j]=c
    
    return splits, adj, len(splits1s), len(splits2s)


#================================================================================#  
#Methods for finding the geodesic path
#================================================================================#  

def geodesic(adj,bl1,bl2,neg,todo): #returns the last orthant


    bestp=[None]
    mind=[cone(bl1,bl2)+1]
    countp=[0]
    
    def best_path(edge,past=[]): #finds shortest path through graph
        past.append(edge)
        if not edge.anc.edges:
            countp[0]+=1
            newp=Path(past)
            dist=newp.distance(bl1,bl2)
            if sum(dist)<mind[0]:
                mind[0]=sum(dist)
                bestp[0]=newp
        for e in edge.anc.edges:
            if e.s < past[-1].s:best_path(e,deepcopy(past))
            
    def one_path(edge):
        elist=[edge]
        while edge.anc.edges:
            edge=edge.anc.edges[0]
            elist.append(edge)
        return Path(elist)

    def next_sub(str,i): return [x for (pos,x) in zip(range(len(str)), str) if (2**pos) & i] #enumerate all 2**str subsets by binary numbers

    #set class variables
    Orthant.opt=True
    Orthant.dim1=len(bl1)
    Orthant.dim2=len(bl2)
    Orthant.adj=adj

    length=2**Orthant.dim2  #number of indices generated for the orthants
    Orthant.OrthList=[None]*length

    ###Orthant.IDList=[-1]*length #stores mapping of simple IDs (binary numbers) to complex IDs (binomial coefficients)
    Orthant.IDList={}
    
    start=Orthant(Neg=neg,todoPos=todo)
    startid=start.get_id()
    Orthant.OrthList[startid]=start

    max_i=0 #index of last orthant

    for i in range(startid,length):
        if not Orthant.OrthList[i]:continue #many orthants are not generated bec suborthants of a larger one
        max_i=i
        orth=Orthant.OrthList[i]

        actdim=len(orth.get_todo())
    
        for j in range(1,2**actdim): #0. Subset is empty
            switch=next_sub(orth.get_todo(),j) 
            newedge=orth.clone(switch) #switch is the initial R, the complete edge is generated in the constructor

            oid=newedge.get_id()
            oldorth=Orthant.OrthList[oid]
            if oldorth and oldorth.edges and (newedge in oldorth.edges):
                del newedge
                continue  #path already computed

            if not oldorth: #generate orthant where the edge points to
                neworth=newedge.create_orthant()
                succ=newedge.compute_s(bl1,bl2,neworth)
                if succ: Orthant.OrthList[oid]=neworth #may be not successfull, because transition times have to fullfill several constraints
                else:
                    del newedge
                    del neworth
                
            else:
                succ=newedge.compute_s(bl1,bl2,oldorth)
                if not succ:continue
                oldorth.edges.append(newedge)
                
                oldedge=oldorth.edges[0]
                maxsk=max(newedge.s,oldedge.s)
                distdiff=sum(one_path(newedge).distance(bl1,bl2))-sum(one_path(oldedge).distance(bl1,bl2))
                if distdiff<0: oldorth.edges=[newedge]
                    
    lastedge=IEdge(0,0,0,0,anc=Orthant.OrthList[max_i],s=1)
    best_path(lastedge)

    return bestp[0],countp[0]

    


#================================================================================#
#Methods for computing all the distances
#================================================================================#  

def cone(diff1,diff2,shared1=[],shared2=[]):
    sharednorm=snorm([shared1[i]-shared2[i] for i in range(0,len(shared1))])
    return math.sqrt((norm(diff1)+norm(diff2))**2+sharednorm)


def distance(tree1,tree2,outfile=None):
    

    def permutation_indices(data): 
        return sorted(range(len(data)), key = data.__getitem__)

    def branch_score(diff1,diff2,shared1=[],shared2=[]):
        shared=[shared1[i]-shared2[i] for i in range(0,len(shared1))]
        return norm(diff1+diff2+shared)

    def combine(diff,shared,splits1,bl1,dstart,dend): #combine branch length lists so that indicees correspond to splits
        shared_branch=[bl1[splits1.index(s)] for s in shared]
        diff_branch=[bl1[splits1.index(diff[i])] for i in range (dstart,dend)]
        return diff_branch,shared_branch
        
    def create_shared_equ(branch1,branch2):
        equs=[]
        for i in range(0,len(branch1)):
            equs.append([0,branch2[i]-branch1[i],branch1[i]])
        return equs
    
    def inverse(mat):return [[mat[i][j] for i in range(0,len(mat))] for j in range(0,len(mat[0]))]

    #================================================================================#  
        
    split_decomp=[[list(get_splits(tree1,True)),list(get_splits(tree2,True))]]
    header=False
    spp=split_decomp[0][0][2]
        
    all_equs=[]
    dec=1 #number of decomposition

    s_dic={} #saves geodesic path in transition times and Left and Right sets
    bl_dic1={} #saves branch length of splits in case of multiple decompositions
    bl_dic2={}
    len_list=["cone","bs","geod","coneall","bsall","geodall"]
    split_list=["shared","diff1","diff2","dim"]
    len_dic={}.fromkeys(len_list) #save length and splits to compute overall numbers in the end
    for l in len_dic.keys():
        len_dic[l]=[0]
    split_dic={}.fromkeys(split_list,0)
    graph=1
    compl_time=0

    for t1,t2 in split_decomp:
        splits1,bl1,spp1=t1
        splits2,bl2,spp2=t2
        new_s_dic={}
        
        # =============================================================================== #
        # creates set S (diff_splits), compatibility matrix (adj), C (shared_splits) and corresponding numbers
 
        diff_splits,adj,dim1,dim2=get_split_representation(splits1,splits2)
    
        shared_splits=list(set(splits1).intersection(set(splits2)))
        shared_splits.sort()
        c=len(shared_splits)
        overlap=set(splits1).intersection(splits2)

        # =============================================================================== #
        # combine branch length lists so that indices correspond to splits
        
        branch1_diff,branch1_shared=combine(diff_splits,shared_splits,splits1,bl1,0,dim1)
        branch2_diff,branch2_shared=combine(diff_splits,shared_splits,splits2,bl2,dim1,dim1+dim2)
    
        # =============================================================================== #
        # output

        split_dic["shared"]+=len(overlap)
        split_dic["diff1"]+=len(branch1_diff)
        split_dic["diff2"]+=len(branch2_diff)
                
        # =============================================================================== #
        # fill bl_dics
        if header:
            for i in range(0,dim1):
                bl_dic1["%u/%u"%(dec,i+1)]=branch1_diff[i]
                bl_dic2["%u/%u"%(dec,i+1)]=0
            for i in range(0,dim2):
                bl_dic1["%u/%u"%(dec,dim1+i+1)]=0
                bl_dic2["%u/%u"%(dec,dim1+i+1)]=branch2_diff[i]
            for i in range (0,c):
                bl_dic1["%u/%u"%(dec,dim1+dim2+i+1)]=branch1_shared[i]
                bl_dic2["%u/%u"%(dec,dim1+dim2+i+1)]=branch2_shared[i]
        else:
            for i in range(0,dim1):
                bl_dic1[str(i+1)]=branch1_diff[i]
                bl_dic2[str(i+1)]=0
            for i in range(0,dim2):
                bl_dic1[str(dim1+i+1)]=0
                bl_dic2[str(dim1+i+1)]=branch2_diff[i]
            for i in range (0,c):
                bl_dic1[str(dim1+dim2+i+1)]=branch1_shared[i]
                bl_dic2[str(dim1+dim2+i+1)]=branch2_shared[i]
            
                
        # =============================================================================== #
        # there are different splits
    
        if dim1 or dim2:
            # =============================================================================== #
            # compute bounds          
            
            bs=branch_score(branch1_diff,branch2_diff) #lower bound
            ub=cone(branch1_diff,branch2_diff) #upper bound
            len_dic["cone"].append(ub)
            len_dic["bs"].append(bs)
            
            
          
            # =============================================================================== #
            # find splits that are compatible to all others

            l_ind=[]
            #There may be some splits in the first tree that are compatible to all splits in the second tree -> they have to end at 1
            for i in range(0,dim1):
                if adj[i]==[True]*dim2:l_ind.append(i)
            comp_equ_l=[[1,-branch1_diff[i],branch1_diff[i]] for i in l_ind]
            neg=set(range(0,dim1)).difference(set(l_ind))  #neg is the starting neg for the geodesic algorithm, in case of full compatibilities, some are excluded

            r_ind=[]
            #vice versa
            for j in range(0,dim2):
                comp=True
                for i in range(0,dim1):
                    if not adj[i][j]:
                        comp=False
                        break
                if comp:r_ind.append(j)
            comp_equ_r=[[0,branch2_diff[i],0] for i in r_ind]
            todo=set(range(0,dim2)).difference(set(r_ind))        
    
            # =============================================================================== #
            # geodesic distance algorithm if still something todo

            t0=time.clock()
            if len(todo)>0: #it may be that all were compatible bec of polytomies

                # =============================================================================== #
                # todo should be the smaller set since algorithm is exponential in len(todo)
                swap=(len(neg)<len(todo))
                if swap:
                    adj=inverse(adj)
                    branch1_diff,branch2_diff=branch2_diff,branch1_diff
                    neg,todo=todo,neg
                
                split_dic["dim"]+=len(todo)

                try:
                    path,count=geodesic(adj,branch1_diff,branch2_diff,neg,todo)
                except OverflowError:
                    print "To many splits in actual decomposition to compute the geodesic distance exactly:", Orthant.dim2
                    print "Start the computation with option '-a' to compute the approximations"
                    continue
                
                t1=time.clock()
           
                equ_l,equ_r= path.equations(branch1_diff,branch2_diff) #the assignment to first and second tree is not correct, but changes nothing for distance computation!!!
                if swap:
                    branch1_diff,branch2_diff=branch2_diff,branch1_diff
                    path.inverse()
        
            else:
                equ_l=[]
                equ_r=[]
                count=1
                path=None
                t1=t0
                

            # =============================================================================== #
            # Get the actual L, R and s and save in s_dic
            
            if path:
                L = path.get_L();L=[[i+1 for i in l] for l in L]
                R = path.get_R();R=[[i+dim1+1 for i in r] for r in R]
                s = path.get_s()
            else:
                L=[[],[]]
                R=[[],[]]
                s=[0,1]
            if r_ind: R[0]+=[i+dim1+1 for i in r_ind]
            if l_ind: L[-1]+=[i+1 for i in l_ind]    
            if header:
                L=[["%u/%u"%(dec,i) for i in l] for l in L]
                R=[["%u/%u"%(dec,i) for i in r] for r in R]

            for i in range(0,len(s)):
               new_s_dic[s[i]]=[L[i],R[i]]


            # =============================================================================== #
            # Compute the distances for different splits and output them

            equ_l+=comp_equ_l
            equ_r+=comp_equ_r

            dist=dist_for_geod(equ_r,equ_l)
            timeg=t1-t0
            

                
            len_dic["geod"].append(sum(dist))
            
        # =================================================================================== #
        # There are no different splits
        
        else:
            timeg=0
            equ_r=[];equ_l=[];count=1 #count is number of geodesic paths
            new_s_dic[0]=[[],[]]
            new_s_dic[1]=[[],[]]
            

        # =============================================================================== #
        # Compute the distances for all splits and output them
        
        bs_compl=branch_score(branch1_diff,branch2_diff,branch1_shared,branch2_shared)
        cone_compl=cone(branch1_diff,branch2_diff,branch1_shared,branch2_shared)
        len_dic["coneall"].append(cone_compl)
        len_dic["bsall"].append(bs_compl)

        
        dist_compl=dist_for_geod(equ_r+create_shared_equ(branch1_shared,branch2_shared),equ_l)
        
        len_dic["geodall"].append(sum(dist_compl))


        for s in new_s_dic.keys():
            if s_dic.has_key(s):
                s_dic[s][0]+=new_s_dic[s][0]
                s_dic[s][1]+=new_s_dic[s][1]
            else: s_dic[s]=new_s_dic[s]
  
        compl_time+=timeg
        dec+=1

    #===============================================================================#
    #output in case of decompositions

    if header:
        for i in len_dic.keys():
            len_dic[i]=norm(len_dic[i])
        
        
        
    return len_dic["geodall"][1]

#===============================================================================#


def main(tree1, tree2):
    
    return distance(tree1[:-1],tree2[:-1],None)
    