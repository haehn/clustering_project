from itertools import permutations
import random, string

c = [1,2,2,3,2,3,1,1,1]
t = [2,2,2,3,3,3,1,1,1]


def sc(c, t):
    def perms(t):
        if type(t)==type([1]):
            t = ''.join([str(x) for x in t])
        a = string.lowercase
        s = sorted(list(set([str(item) for item in t])))
        d = dict(zip(s,a))
        per = [list(x) for x in permutations(s)]
        sub = ''.join([str(item) for item in t])
        for x in s:
            sub=sub.replace(x,d[x])
        l=[]
        for poss in per:
            d2 = dict(zip(a,poss))
            new = sub
            for x in poss:
                new = new.replace(d[x],d2[d[x]])
            l.append(new)
        l= [ list(x) for x in l]
        return l
    
    def comp(l,l2):
        score=0
        for (x,y) in zip(l,l2):
            if x != y : score+=1
        return score

    if type(c)==type([1]):
        c = ''.join([str(x) for x in c])
    return min(map(comp, [c]*len(perms(t)), perms(t)))

print sc(c,t)
