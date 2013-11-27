# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# *Keywords* dataset

# <markdowncell>

# [Luís F. Simões](mailto:luis.simoes@vu.nl)<br>
# 2013-11-26<div style="float: right">`Notebooks:` [index &uarr;](https://github.com/lfsimoes/VU/tree/master/2013__Collective_Intelligence)</div><br><br>
# 
# *******

# <markdowncell>

# As you can see on the [Building the "evolution" research papers dataset](http://nbviewer.ipython.org/urls/raw.github.com/lfsimoes/VU/master/2013__Collective_Intelligence/11__Entrez__building_dataset.ipynb#EFetch:-Downloading-full-records-from-Entrez) notebook, one of the fields retrieved by `Entrez.efetch()` is a `KeywordList`. Keywords are sometimes provided by authors, next to papers' abstracts, as a way of categorizing the research work presented in the paper. Here we will build a *Keywords* dataset, mapping each PubMed ID in our core dataset, to its list of keywords (if available).

# <codecell>

from Bio import Entrez

# NCBI requires you to set your email address to make use of NCBI's E-utilities
Entrez.email = "Your.Name.Here@example.org"

# <codecell>

import os, cPickle, bz2
from collections import *

# data format used in Python object serialization
pickle_protocol = 2

# <codecell>

search_term = 'evolution'

# <codecell>

Ids_file = search_term + '__Ids.pkl.bz2'

Ids = cPickle.load( bz2.BZ2File( Ids_file, 'rb' ) )

# <codecell>

def conv_str( s ):
    return unicode(s) if isinstance(s, unicode) else str(s)

# <markdowncell>

# *******

# <codecell>

Keywords_file = search_term + '__Keywords.pkl.bz2'

# <codecell>

import httplib, urllib2
from collections import deque

if os.path.exists( Keywords_file ):
    Keywords = cPickle.load( bz2.BZ2File( Keywords_file, 'rb' ) )
else:
    # `Keywords` will be incrementally assembled, by performing multiple queries,
    # each returning at most `retrieve_per_query` entries.
    Keywords = deque()
    retrieve_per_query = 200
    
    print 'Fetching Keywords of results: ',
    for start in xrange( 0, len(Ids), retrieve_per_query ):
        print start,
        
        # build comma separated string with the ids at indexes [start, start+retrieve_per_query)
        query_ids = ','.join( [ str(id) for id in Ids[ start : start+retrieve_per_query ] ] )
        
        # issue requests to the server, until we get the full amount of data we expect
        while True:
            try:
                s = Entrez.read( Entrez.efetch(db="pubmed", id=query_ids, retmode="xml" ) )
            except (httplib.IncompleteRead, httplib.BadStatusLine, urllib2.HTTPError):
                print 'r',
                continue
            break
        
        for p in s:
            key  = 'MedlineCitation' if 'MedlineCitation' in p else 'BookDocument'
            pmid = p[key]['PMID']
            kws  = p[key]['KeywordList']
            if kws != []:
                kws  = [ conv_str(kw).lower() for kw in kws[0] ]
            
            Keywords.append( (int(pmid), kws) )
    
    # Save Keywords, as a dictionary indexed by Ids
    Keywords = dict( Keywords )
    
    cPickle.dump( Keywords, bz2.BZ2File( Keywords_file, 'wb' ), protocol=pickle_protocol )

# <markdowncell>

# Building the reverse mapping:

# <codecell>

from collections import defaultdict

papers_with_keyword = defaultdict(set)

for pid,ks in Keywords.iteritems():
    for k in ks:
        papers_with_keyword[ k ].add( pid )

# <markdowncell>

# Basic dataset statistics:

# <codecell>

with_kws = sum(1 for k in Keywords.itervalues() if k!=[])
print '%d papers (%.2f %%) have keywords data' % (with_kws, 100.*with_kws/len(Keywords))

# <codecell>

print '%d keywords used in total' % len(papers_with_keyword)

# <markdowncell>

# Finding the 50 most used keywords:

# <codecell>

uses_per_keyword = sorted( [ (k,len(v)) for k,v in papers_with_keyword.iteritems()], key=lambda i:i[1] )
uses_per_keyword[-50:]

# <markdowncell>

# Finding the paper tagged with the greatest number of keywords:

# <codecell>

max( Keywords.iteritems(), key=lambda i:len(i[1]) )

# <markdowncell>

# The PubMed ID 24012599 corresponds to the paper [Polishing the craft of genetic diversity creation in directed evolution](http://www.ncbi.nlm.nih.gov/pubmed/?term=24012599).

# <markdowncell>

# *******

# <codecell>

import matplotlib.pyplot as plt

# show plots inline within the notebook
%matplotlib inline
# set plots' resolution
plt.rcParams['savefig.dpi'] = 100

# <markdowncell>

# How many papers use a given number of keywords?
# 
# Only papers containing at least one keyword are considered in the plot below.

# <codecell>

plt.hist( [ len(ks) for ks in Keywords.itervalues() if len(ks)>0 ], bins=range(60), align='left' )
plt.xlabel( 'Number of keywords in a certain paper' )
plt.ylabel( 'Number of papers' )
#plt.yscale('log')
plt.xlim(0,50);

# <markdowncell>

# How many keywords are used in a given number of papers?

# <codecell>

plt.hist( [kwus[1] for kwus in uses_per_keyword], bins=range(25), histtype='step', align='left' )
plt.xlabel( 'Number of papers using a certain keyword' )
plt.ylabel( 'Number of keywords' )
plt.yscale('log')
plt.xlim(0,20);

# <markdowncell>

# We see above, for instance, that there are ~1000 keywords which are used in only 3 papers.

# <markdowncell>

# *******

