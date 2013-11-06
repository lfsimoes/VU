# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Network analysis

# <markdowncell>

# [Luís F. Simões](mailto:luis.simoes@vu.nl)<br>
# 2013-11-06<div style="float: right">`Notebooks:` [&larr; previous](http://nbviewer.ipython.org/urls/raw.github.com/lfsimoes/VU/master/2013__Collective_Intelligence/12__inspecting_the_data.ipynb) &bull; [index &uarr;](https://github.com/lfsimoes/VU/tree/master/2013__Collective_Intelligence)</div><br><br>
# 
# *******

# <markdowncell>

# This notebook exemplifies ways of using those parts of the dataset that are better seen as representing graphs. Namely, the co-authorship, and citations networks.

# <headingcell level=2>

# Loading modules & datasets

# <codecell>

import cPickle, bz2
from collections import *


import numpy as np
import matplotlib.pyplot as plt

# show plots inline within the notebook
%matplotlib inline
# set plots' resolution
plt.rcParams['savefig.dpi'] = 100

from IPython.display import display, HTML

# <codecell>

search_term = 'evolution'

Ids_file       = search_term + '__Ids.pkl.bz2'
Summaries_file = search_term + '__Summaries.pkl.bz2'
Citations_file = search_term + '__Citations.pkl.bz2'

# <codecell>

Ids = cPickle.load( bz2.BZ2File( Ids_file, 'rb' ) )

# <codecell>

Summaries = cPickle.load( bz2.BZ2File( Summaries_file, 'rb' ) )

paper = namedtuple( 'paper', ['title', 'authors', 'year', 'doi'] )

for (id, paper_info) in Summaries.iteritems():
    Summaries[id] = paper( *paper_info )

# <codecell>

Citations = cPickle.load( bz2.BZ2File( Citations_file, 'rb' ) )

# <headingcell level=2>

# Auxiliary functions

# <markdowncell>

# The function below will allow us to display paper summaries in an IPython-aware way.
# 
# You can find examples [here](http://nbviewer.ipython.org/urls/raw.github.com/ipython/ipython/1.x/examples/notebooks/Part%205%20-%20Rich%20Display%20System.ipynb) of how to use IPython's Rich Display System, or the whole set of examples [here](https://github.com/ipython/ipython/tree/1.x/examples/notebooks).

# <codecell>

def display_summary( id, extra_text='' ):
    """
    Function for printing a paper's summary through IPython's Rich Display System.
    Trims long titles or author lists, and links to the paper's  DOI (when available).
    """
    s = Summaries[ id ]
    
    title = ( s.title if s.title[-1]!='.' else s.title[:-1] )
    title = title[:150].rstrip() + ('' if len(title)<=150 else '...')
    if s.doi!='':
        title = '<a href=http://dx.doi.org/%s>%s</a>' % (s.doi, title)
    
    authors = ', '.join( s.authors[:5] ) + ('' if len(s.authors)<=5 else ', ...')
    
    lines = [
        title,
        authors,
        str(s.year),
        '<small>id: %d%s</small>' % (id, extra_text)
        ]
    
    display( HTML( '<blockquote>%s</blockquote>' % '<br>'.join(lines) ) )

# <codecell>

# usage example
display_summary( 23144668 )

# <headingcell level=2>

# Co-authorship network

# <markdowncell>

# `Summaries` maps paper *ids* to paper *summaries*. Let us now create here mappings by different criteria.
# 
# We'll start by building a mapping from *authors*, to *ids* of papers they authored:

# <codecell>

papers_of_author = defaultdict(list)

for id,p in Summaries.iteritems():
    for a in p.authors:
        papers_of_author[ a ].append( id )

# <codecell>

papers_of_author['Eiben AE']

# <codecell>

for id in papers_of_author['Eiben AE']:
    display_summary( id )

# <markdowncell>

# We now build a co-authorship network, a graph linking *authors*, to the set of *co-authors* they have published with.
# 
# To obtain the union of author lists, we'll be using Python's [sets](http://docs.python.org/2/library/stdtypes.html#set-types-set-frozenset).

# <codecell>

coauthors = defaultdict(set)

for p in Summaries.itervalues():
    for a in p.authors:
        coauthors[ a ].update( p.authors )

# the code above results in each author being listed as having co-autored with himself. We now remove such references here
for a,ca in coauthors.iteritems():
    ca.remove( a )

# <codecell>

print ', '.join( coauthors['Eiben AE'] )

# <headingcell level=4>

# [Degree distribution](https://en.wikipedia.org/wiki/Degree_distribution)

# <markdowncell>

# With this data in hand, we can plot the distribution showing the number of collaborators a scientist has published with in its full publication record:

# <codecell>

plt.hist( x=[ len(ca) for ca in coauthors.itervalues() ], bins=range(55), histtype='bar', align='left', normed=True )
plt.xlabel('number of collaborators')
plt.ylabel('fraction of scientists')
plt.xlim(0,50);

# <headingcell level=2>

# Citations network

# <markdowncell>

# We'll start by expanding the `Citations` dataset into two mappings: 
# 
# * `papers_citing[ id ]`: papers citing a given paper;
# * `cited_by[ id ]`: papers cited by a given paper (in other words, its list of references).
# 
# If we see the Citations dataset as a directed graph where papers are nodes, and citations links between then, then `papers_citing` gives you the list of a node's incoming links, whereas `cited_by` gives you the list of its outgoing links.
# 
# The dataset was assembled by querying for papers citing a given paper. As a result, the data mapped to in `cited_by` (its values) is necessarily limited to ids of papers that are part of the dataset (they are members of `Ids`). Two papers that are never cited, but do themselves cite many of the same papers (papers that are not members of `Ids`), are probably similar/relevant to each other. Unfortunately, that kind of data is not currently available.

# <codecell>

papers_citing = Citations  # no changes needed, this is what we are storing already in the Citations dataset

cited_by = defaultdict(list)

for ref, papers_citing_ref in papers_citing.iteritems():
    for id in papers_citing_ref:
        cited_by[ id ].append( ref )

# <markdowncell>

# Back in the first notebook we used as an example the paper with PubMed ID 20949101 ("Critical dynamics in the evolution of stochastic strategies for the iterated prisoner's dilemma"). We can now use the `cited_by` mapping to retrieve what we know of its list of references.
# 
# As mentioned above, because our server queries asked for papers citing a given paper (and not papers a paper cites), the papers we get through `cited_by` are then necessarily all members of our datasets, and we can therefore find them in `Summaries`.
# 
# You can match the list below with [the paper's References section](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2951343/#__ref-listid953859title). We have identified 22 of its references, from a total of 52.

# <codecell>

paper_id = 20949101
refs = { id : Summaries[id].title for id in cited_by[20949101] }
print len(refs), 'references identified for the paper with id', paper_id
refs

# <markdowncell>

# If we lookup the same paper in `papers_citing`, we now see that one of the papers citing it is itself part of our dataset, but the other one (with PubMed ID 23903782) is instead a citation coming from outside the dataset. Back in the first notebook, we saw that this id corresponds to the paper 'Evolutionary instability of zero-determinant strategies demonstrates that winning is not everything' (quite an interesting paper actually). But, being an "outsider", its information is not in `Summaries`. Of course, you can always issue a new Entrez query to get its details if you want to.

# <codecell>

{ id : Summaries.get(id,['??'])[0]  for id in papers_citing[paper_id] }

# <markdowncell>

# Of course, though no direct query was ever issued for this paper id 23903782, its repeated occurrence in other papers' citation lists does allow us to reconstruct a good portion of its references. Below is the list of papers in our dataset cited by that paper (10 references, out of the paper's total of 32, as you can see [here](http://www.nature.com/ncomms/2013/130801/ncomms3193/full/ncomms3193.html#references)):

# <codecell>

paper_id2 = 23903782
refs2 = { id : Summaries[id].title for id in cited_by[paper_id2] }
print len(refs2), 'references identified for the paper with id', paper_id2
refs2

# <markdowncell>

# We can easily combine boh mappings (incoming and outgoing links) to answer such questions as: "which papers are the most cited by those that cite a given paper?".
# 
# Going back to our example paper ("Critical dynamics in the evolution of stochastic strategies for the iterated prisoner's dilemma"), we find two references that are commonly cited by the two papers citing it. If you look above to the list of papers our paper cites, you find that actually, the paper itself also cites these same two references, indicating that they are somehow particularly relevant to this line of research.

# <codecell>

# Counter built over list of papers cited by papers that cite a given id
cocited = Counter([
    ref
    for citers in papers_citing[ paper_id ]
    for ref in cited_by[ citers ]
    if ref != paper_id
    ])

# discard papers cited only once
cocited = { id : nr_cocits for (id, nr_cocits) in cocited.iteritems() if nr_cocits > 1 }


for (id, times_co_cited) in sorted( cocited.items(), key=lambda i:i[1], reverse=True ):
    display_summary( id, ', nr. co-citations: %d' % times_co_cited )

# <markdowncell>

# A final implementation note regarding `cited_by`: the way it was built, only papers cited at least once appear as keys in it. With `cited_by` being a [defaultdict](http://docs.python.org/2/library/collections.html#collections.defaultdict), attempts to access a missing element will silently cause it to be added (with an empty list as vale, in this case), instead of raising an exception. Just take note that if you iterate over keys in `cited_by`, you are not iterating over all of the graph's nodes.
# 
# As an example, here's an isolated node in the graph, an old paper for which we have no citation data (neither to nor from).

# <codecell>

display_summary( Ids[-1] )

# <codecell>

papers_citing[ 21373130 ]

# <codecell>

21373130 in cited_by

# <codecell>

cited_by[ 21373130 ]

# <codecell>

21373130 in cited_by

# <codecell>

# just to revert to the way things were
del cited_by[ 21373130 ]

# <markdowncell>

# Should you want to "lock" `cited_by`, to make it behave as a regular Python dictionary, just convert it like this:

# <codecell>

cited_by = dict(cited_by)

# <headingcell level=3>

# Dataset statistics

# <markdowncell>

# Now that we have a better understanding about the data we're dealing with, let us obtain some basic statistics about our graph.

# <codecell>

print 'Number of core ids %d (100.00 %%)' % len(Ids)

with_cit = [ id for id in Ids if papers_citing[id]!=[] ]
print 'Number of papers cited at least once: %d (%.2f %%)' % (len(with_cit), 100.*len(with_cit)/len(Ids))

isolated = set( id for id in Ids if papers_citing[id]==[] and id not in cited_by )
print 'Number of isolated nodes: %d (%.2f %%)\n\t'   \
      '(papers that are not cited by any others, nor do themselves cite any in the dataset)'% (
    len(isolated), 100.*len(isolated)/len(Ids) )

noCit_withRefs = [ id for id in Ids if papers_citing[id]==[] and id in cited_by ]
print 'Number of dataset ids with no citations, but known references: %d (%.2f %%)' % (
    len(noCit_withRefs), 100.*len(noCit_withRefs)/len(Ids))

print '(percentages calculated with respect to just the core ids (members of `Ids`) -- exclude outsider ids)\n'

# <codecell>

Ids_set    = set( Ids )
citing_Ids = set( cited_by.keys() ) # == set( c for citing in papers_citing.itervalues() for c in citing )

outsiders = citing_Ids - Ids_set    # set difference: removes from `citing_Ids` all the ids that occur in `Ids_set`
nodes     = citing_Ids | Ids_set - isolated     # set union, followed by set difference

print 'Number of (non-isolated) nodes in the graph: %d\n\t(papers with at least 1 known citation, or 1 known reference)' % len(nodes)
print len( citing_Ids ), 'distinct ids are citing papers in our dataset.'
print 'Of those, %d (%.2f %%) are ids from outside the dataset.\n' % ( len(outsiders), 100.*len(outsiders)/len(citing_Ids) )

# <codecell>

all_cits      = [ c for citing in papers_citing.itervalues() for c in citing ]
outsider_cits = [ c for citing in papers_citing.itervalues() for c in citing if c in outsiders ]

print 'Number of links (citations) in the graph:', len(all_cits)
print 'A total of %d citations are logged in the dataset.' % len(all_cits)
print 'Citations by ids from outside the dataset comprise %d (%.2f %%) of that total.\n' % (
    len(outsider_cits),
    100.*len(outsider_cits)/len(all_cits) )

# <headingcell level=4>

# Most cited papers

# <markdowncell>

# Let us now find which are the 10 most cited papers in our dataset.

# <codecell>

nr_cits_per_paper = [ (id, len(cits)) for (id,cits) in papers_citing.iteritems() ]

for (id, cits) in sorted( nr_cits_per_paper, key=lambda i:i[1], reverse=True )[:10]:
    display_summary( id, ', nr. citations: <u>%d</u>' % cits )

# <markdowncell>

# The most cited paper in our dataset "Initial sequencing and analysis of the human genome", has 3258 citations, which is ~20% of the citations Google Scholar recognizes for it (see [here](http://scholar.google.com/scholar?q=10.1038%2F35057062)).

# <headingcell level=4>

# [Degree distribution](https://en.wikipedia.org/wiki/Degree_distribution): Citations (indegree distribution)

# <codecell>

nr_cits = [ len(cits) for cits in papers_citing.itervalues() if cits!=[] ]

# <codecell>

plt.hist( nr_cits, bins=range(103), histtype='bar', align='left' )
plt.xlabel('number of citations')
plt.ylabel('number of papers')
plt.xlim(0,100);

# <codecell>

plt.hist( nr_cits, bins=range(103), histtype='step', align='left', normed=True )
plt.xlabel('number of citations')
plt.ylabel('fraction of papers')
plt.yscale('log')
#plt.xscale('log')
plt.xlim(0,100);

# <headingcell level=4>

# [Degree distribution](https://en.wikipedia.org/wiki/Degree_distribution): References (outdegree distribution)

# <codecell>

nr_refs = [ len(refs) for refs in cited_by.itervalues() ]

# <codecell>

plt.hist( nr_refs, bins=range(23), histtype='bar', align='left' )
plt.xlabel('number of references')
plt.ylabel('number of papers')
plt.xlim(0,20);

# <codecell>

plt.hist( nr_refs, bins=range(23), histtype='step', align='left', normed=True )
plt.xlabel('number of references')
plt.ylabel('fraction of papers')
plt.yscale('log')
#plt.xscale('log')
plt.xlim(0,20);

