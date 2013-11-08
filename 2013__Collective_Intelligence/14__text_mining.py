# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Text mining

# <markdowncell>

# [Luís F. Simões](mailto:luis.simoes@vu.nl)<br>
# 2013-11-08<div style="float: right">`Notebooks:` [previous &larr;](http://nbviewer.ipython.org/urls/raw.github.com/lfsimoes/VU/master/2013__Collective_Intelligence/13__network_analysis.ipynbb) &bull; [index &uarr;](https://github.com/lfsimoes/VU/tree/master/2013__Collective_Intelligence)</div><br><br>
# 
# *******

# <headingcell level=2>

# Loading modules & datasets

# <codecell>

import cPickle, bz2
from collections import *

# <codecell>

search_term = 'evolution'

Summaries_file = search_term + '__Summaries.pkl.bz2'
Abstracts_file = search_term + '__Abstracts.pkl.bz2'

# <codecell>

Summaries = cPickle.load( bz2.BZ2File( Summaries_file, 'rb' ) )

paper = namedtuple( 'paper', ['title', 'authors', 'year', 'doi'] )

for (id, paper_info) in Summaries.iteritems():
    Summaries[id] = paper( *paper_info )

# <codecell>

Abstracts = cPickle.load( bz2.BZ2File( Abstracts_file, 'rb' ) )

# <headingcell level=2>

# Abstracts dataset

# <codecell>

with_abstr = sum( 1 for a in Abstracts.itervalues() if a!='' )
print 'Abstracts available for %d (%.2f %%) of the papers in the dataset.' % (with_abstr, 100.*with_abstr/len(Abstracts))

# <headingcell level=3>

# Word frequencies

# <markdowncell>

# The function below will allows us to identify words occurring in a string:

# <codecell>

def words( text ):
    return Counter([
        word.replace(',','').replace('.','').lower()
        for word in text.split(' ')
        if word != ''     # the split by spaces generates empty strings when consecutive spaces occur in the title; this discards them
        ])

# <markdowncell>

# Here we apply it to two abstracts:

# <codecell>

Abstracts[ 7466396 ]

# <codecell>

w1 = words( Abstracts[ 7466396 ] )
w1

# <codecell>

w2 = words( Abstracts[ 8316296 ] )
w2

# <markdowncell>

# We can intersect the sets of words occurring in each abstract to obtain those that appear in both.

# <codecell>

shared = sorted( set( w1.keys() ) & set( w2.keys() ), key=lambda i:w1[i]+w2[i], reverse=True )
shared

# <markdowncell>

# And we can now represent each abstract just by this set of selected words:

# <codecell>

[
    [ bag.get(w,0) for w in shared ]
    for bag in [w1,w2]
    ]

