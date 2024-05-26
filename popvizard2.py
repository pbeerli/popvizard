#!/usr/bin/env python
############################################################################
# Population simulations using standard models with exponential
# growth/shrinkage
# Python code creates a single PDF page with a population through time plot
# and a genealogy of a sample.
# Copyright 2011 Peter Beerli
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http:
############################################################################

import sys
import scipy as sc 
import numpy as np
import time
import ast
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


############################################################################
#DEFAULTS
############################################################################
GENERATIONS  = 50
NE           =  30
GROWTH       = 0.00
BOTTLENECK   = "{0:30}"
# samples taken from population
SAMPLESIZE   = 10
#
SEED         = None
MODEL        = 'WrightFisher' 
VAR          = 1.0
MSIZE        = 3.0
SCOLOR       = 'b'
FILENAME     = 'wf.pdf' 

############################################################################

def intersection(c1,c2):
    return c2 in c1


#############################################################################


def mymodel(start, ne, oldne, std, RS, mymodel='Moran'):
    if mymodel=='Moran':
        diff = ne - oldne
        #print "oldne", oldne
        #print "ne", ne
        #print "len(start)",len(start)
        #print "diff",diff
        a=list(range(0,ne))
        if diff>0:
            a.pop(RS.randint(0,len(a)))
            for i in range(diff+1):
                a.append(RS.randint(0,ne))
        elif diff == 0:
            a.pop(RS.randint(0,len(a)))
            a.append(RS.randint(0,ne))   
        else:
            for i in range(-diff):
                a.append(RS.randint(0,ne))
            a.pop(RS.randint(0,len(a)))
            a.append(RS.randint(0,ne))
    else:
        if mymodel=='WrightFisher':
            a = RS.randint(0,ne,oldne)
        else: # mymodel=='Canning':
            alpha = 1.0 / (std * std)
            beta = 1./alpha
            offsprings=[]
            for x in range(0,ne):
                limi = 1+int(RS.gamma(alpha,beta,1))
                #print "limit",limi
                offsprings.append([x for i in range(0,limi)])
            offsprings = reduce(lambda x,y: x+y,offsprings)
            while len(offsprings)<oldne:
                offsprings.append(offsprings[RS.randint(0,len(offsprings))]) 
                    
            #print "offsprings", len(offsprings), offsprings 
            RS.shuffle(offsprings) 
            a = offsprings[:oldne]
            #print "len(a)",len(a) 
    a.sort()
    return a


def popsim(generations=50,ne=30,growth=0,bottleneck={},samples=[], seed=None, model='Moran', std=1.0, msize=3.0, scolor='b',filename=FILENAME,mydpi=0,myratio=[0,0],mysize=[0,0]):
    RS = np.random.mtrand.RandomState()
    if seed == None:
        seed = int(time.time())
    RS.seed(int(seed))
    samplesize = len(samples)
    s = samples
    ne0 = ne
    #print ne, ne0, bottleneck
    start = np.array(list(range(0,ne0)))
    #print "ge",generations
    #print "gr", growth
    #print "ne0", ne0
    ne1 = int(np.exp(-growth * float(generations)) * float(ne0))
    try:
        if bottleneck[0] < ne0:
            ne0 = bottleneck[0]
    except:
        pass
    for bi in bottleneck:
        if ne1 < bottleneck[bi]:
            ne1 = bottleneck[bi]

    if ne1> ne0:
        nemax=ne1
    else:
        nemax=ne0
    #print "ne1",ne1
    #print np.ones(ne1)
    #print start
    end = (np.array([np.ones(ne1)*(generations),np.array(list(range(0,ne1)))])).T 
    x = []
    cgray=[]
    ccolor=[]
    cpointcolor=[]
    cpointgray=[]
    
    #print model
    #print ne
    #print start
    ne=ne0
    for i in range(0,generations):
        oldne = ne
        start = np.array(list(range(0,oldne)))
        #print oldne,i
        if i in bottleneck:
            ne0 = bottleneck[i]
#        if len(bottleneck)>1: 
#            start = np.array(range(0,bottleneck[0]))
        ne = int(np.exp(-growth * i) * ne0)
        #print "@@@", ne, oldne
        a = mymodel(start,ne,oldne, std, RS,model)
        #print "a", len(a),a 
        #print "start",len(start),start
        b = [[[i , start[j]],[i+1,a[j]]] for j in range(0,oldne)]
        for j in b:
            if intersection(s,j[0][1])==False:
                cpointgray.append(j[0])
                cgray.append(j[0])
                cgray.append(j[1])
                cgray.append([None,None])
            else:
                cpointcolor.append(j[0])
                ccolor.append(j[0])
                ccolor.append(j[1])
                ccolor.append([None,None])
        snew=[]
        for j in b:
            if intersection(s,j[0][1])==True:
                snew.append(j[1][1])
        s = snew[:]
        x.append([cgray,ccolor])
    count = 0
    for e1 in end:
        count += 1
        if count > ne:
            break
        if intersection(s,e1[1])==False:
            cpointgray.append(e1)
        else:
            cpointcolor.append(e1)
    
    msizeS=1.2*msize
    ############################################################################
    # create figure for discrete population plot
    #
    #
    fig = plt.figure()
    #mydpi=200
    #mysize = [8.0,11.5]
    
    
    if mydpi > 0:
        fig.set_dpi(mydpi)
    DPI = fig.get_dpi()
    print("DPI:               ", DPI)

    if mysize[0] > 0 and mysize[1] > 0:
            fig.set_size_inches(mysize[0],mysize[1])
    DefaultSize = fig.get_size_inches()
    print ("Default size:      ", DefaultSize[0],"x", DefaultSize[1], "in*in")
    print ("Image size:         %i x %i pt*pt" %(DPI*DefaultSize[0], DPI*DefaultSize[1]))
    if (myratio[0]==0) and (myratio[1]==0):
        myratio = [nemax/(samplesize+2),2]
    gs = gridspec.GridSpec(2, 1,height_ratios=myratio)
    # set up subplots
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])

    ax1.set_clip_on(False)
    ax1.set_frame_on(True)
    graylines= ax1.plot([x[0] for x in cgray],[y[1] for y in cgray],color='0.75')
    plt.setp(graylines,linewidth=1.0)
    ax1.plot([x[0] for x in cpointgray],[y[1] for y in cpointgray],color='0.75',marker='o',linestyle='None',markersize=msize,markeredgecolor='0.75')
    ax1.plot([x[0] for x in ccolor],[y[1] for y in ccolor],color=scolor)
    ax1.plot([x[0] for x in cpointcolor],[y[1] for y in cpointcolor],color=scolor,marker='o',linestyle='None',markersize=msizeS,markeredgecolor=scolor)
    ax1.set_xlabel('Time [Generations into the Past]')
    if growth==0.0:
        ax1.set_ylabel('%d Individuals' % ne)
    else:
        ax1.set_ylabel('Individuals: $N_e^{(0)}$=%d, $N_e^{(%d)}$=%d' % (ne0,generations,ne))
    ax1.set_xlim(-1.,generations+1)
    ax1.set_ylim(-1.,nemax)
    ###############
    # create tree from the above model
    #       
    x=[]
    y=[]
    for i in ccolor:
        if i == [None, None]:
            y.append(x)
            x=[]
        else:
            x.append(i)
    y.sort()
    lasti=y[0][0][0]
    z=[]
    zall=[]
    for i in y: 
        if i[0][0]==lasti:
            z.append(i)
        else:
            lasti=i[0][0]
            zall.append(z)
            z=[]
            z.append(i)
    tips = zall.pop(0)
    #print "tips", tips
    for t in tips:
        #print "t",t 
        for za in zall:
            for z in za:
                if t[-1] == z[0]:
                    t.append(z[1])
    atips = np.array(tips).T 
    atips = np.delete(atips, np.s_[0], axis=0)
    onegen = np.arange(nemax/(len(sample))/2.,nemax,nemax/(len(sample)))
    tree=[]
    for i in atips[0]:
        #        print onegen
        d = {}
        for a, b in zip(i.tolist(), onegen):
            d.setdefault(a, []).append(b)
        onegen = []
        for key in d:
            #print key, len(d[key]),d[key]
            for di in range(len(d[key])):
                onegen.append(sum(d[key])/len(d[key]))
#print "key", key,i.tolist().index(key),len(d[key]), sum(d[key])/len(d[key]) 
#print
#print "one", onegen
#print
        onegen.sort()
        tree.append(onegen)
    count = 0
    newtree   = [tree[0]]
    times=[count]
    count +=1
    for ti  in range(1,len(tree)):
        t = tree[ti]
        xx = np.array(t) - np.array(newtree[-1])
        if np.all(xx) and len(set(t))>1:
            newtree.append(t)
            times.append(count)
        else:
            newtree.append(tree[ti-1])
            times.append(count)
            newtree.append(t)
            times.append(count)
        #print count, t
        count += 1
    lineages  = np.array(newtree).T 
    ax2.axes.get_yaxis().set_visible(False)
    #
    # plot second figure with tree
    #
    # plots the original crooked tree
    #for t in tips:
    #ax2.plot([x[0] for x in t],[y[1] for y in t],color='r')
    for t in lineages:
        #print t
        #print times
        ax2.plot(times, t,color=scolor)
    plt.savefig(filename,format='pdf')
    return [cpointcolor,cpointgray,cgray, ccolor,[times,lineages]]


#######################################################################
#
# main
#
#######################################################################
if __name__ == '__main__':
    
    
    import argparse

    parser = argparse.ArgumentParser(description=
    '''
    Population simulation for several models including 
    exponentially growing or shrinking
    '''
    )
    parser.add_argument('-m','--model',  default=MODEL, action='store', dest='model', help='set the population model, models other than WRIGHTFISHER, CANNING, and MORAN will fail')
    allowed_models = ('CANNING','WRIGHTFISHER','MORAN')
    
    parser.add_argument('-o','--offspringvar',  default=VAR, action='store', type=float, dest='offspring_variance', help='set the offspring variance, this options has only an effect on the CANNING model, good values are 0.5 or 1.5 etc, values close to 0.0 will result in very long coalescent trees, high values will result in very short coalescent trees')

    parser.add_argument('-n','--Ne',  default=NE, action='store', type=int, dest='ne', help='the effective population size today')

    parser.add_argument('-t','--generations',  default=GENERATIONS, action='store', type=int, dest='generations', help='set the number of generations to plot')
    
    parser.add_argument('-r','--growth',  default=GROWTH, action='store', type=float, dest='growth', help='set the exponential growth rate: -0.01 or 0.01 are good values')

    parser.add_argument('-b','--bottleneck',  default=BOTTLENECK, action='store', type=str, dest='bottleneck', help='set a list of population sizes using python dictionary notation \"{generation1:ne1,generation2:ne2,....}\"')

    parser.add_argument('-s','--samplesize',  default=SAMPLESIZE, action='store', type=int, dest='samplesize', help='sample size')

    parser.add_argument('--seed',  default=SEED, action='store', dest='seed', help='seed for random number generator')

    parser.add_argument('-f', '--filename',  default=FILENAME, action='store', dest='filename', help='filename for output (format is PDF)')

    parser.add_argument('-ms', '--markersize',  default=MSIZE, action='store', dest='markersize', type=float, help='size of marker to plot, for large Ne and large generation times use a smaller value than 3.0')

    parser.add_argument('-c', '--color',  default=SCOLOR, action='store', dest='color', help='color of the coalescent sample tree')
    
    parser.add_argument('-rat', '--ratio',  default=None, action='store', type=float, dest='ratio', help='Ratio between size of population plot and genealogy plot, values larger than 1.0 make the genealogy plot larger than the population plot')

    parser.add_argument('-dpi', '--dpi',  default=None, action='store', type=int, dest='dpi', help='DPI: dots per inch, None means default, perhaps this should be changed for prodcution plots')

    parser.add_argument('-size', '--papersize',  default=[8,11.5], action='store', dest='papersize', help='size of the paper, this needs to be a list of two values, for example [8,11.5]')

    args = parser.parse_args()

    themodel = args.model
    thestd = np.sqrt(args.offspring_variance)
    thegrowth = args.growth
    thebottleneck = ast.literal_eval(args.bottleneck)
    generations  = args.generations
    ne           = args.ne
    # samples taken from population
    samplesize = args.samplesize
    if samplesize>ne:
        samplesize=ne
    population = list(range(0,ne))
    np.random.shuffle(population)
    sample       = population[:samplesize] # samplesize
    sample.sort() # this represents a particular set of individuals
    seed = args.seed
    filename = args.filename
    msize = args.markersize
    scolor = args.color
    dpi = args.dpi   
    ratio = args.ratio
    papersize = args.papersize

    if not(themodel in allowed_models):
        print ("Type", sys.argv[0], "--help for allowed settings")

    print ("Population model:  ", themodel)
    if themodel=='MORAN':
        thestd = np.sqrt(1.0/ne)
    elif themodel=='WRIGHTFISHER':
        thestd= 1.0
    else:
        pass
    print ("Effective size:    ", ne)
    print ("Generations:       ", generations)
    print ("Samplesize:        ", samplesize)
    print ("Offspring variance:", thestd*thestd)
    print ("Growth rate:       ", thegrowth)
    print ("Bottleneck at      ", thebottleneck)
    print ("Random seed:       ", seed)
    print ("Filename:          ", filename)
    print ("Marker size:       ", msize)
    print ("Sample color:      ", scolor)
    print ("Plot ratio:        ", ratio)
    width, height = papersize
    print ("Plot size:         ", width, "x", height)
    if dpi==None:
        thedpi=0
    else:
        thedpi = float(dpi)
    if ratio!=None:
        theratio = [1.0/ratio, 1.0]
    else:
        theratio = [0,0]
    #print "popsim(",generations,ne,thegrowth,thebottleneck, sample,seed, themodel,thestd, msize, scolor,filename,thedpi,theratio, papersize,")"
    cpointcolor,cpointgray,cgray, ccolor,lineages = popsim(generations,ne,thegrowth,thebottleneck, sample,seed, themodel,thestd, msize, scolor,filename,thedpi,theratio, papersize)
    sys.exit(0)
    
