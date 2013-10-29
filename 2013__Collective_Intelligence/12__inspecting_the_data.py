# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Inspecting the dataset

# <markdowncell>

# [Luís F. Simões](mailto:luis.simoes@vu.nl), 2013-10-29<br><br><br>

# <markdowncell>

# This notebook's purpose is to provide a basic illustration of how to handle data in the *"evolution"* dataset.

# <headingcell level=2>

# Loading dataset

# <codecell>

search_term = 'evolution'

Ids_file = search_term + '__Ids.pkl.bz2'
Summaries_file = search_term + '__Summaries.pkl.bz2'

# <codecell>

import cPickle, bz2

Ids = cPickle.load( bz2.BZ2File( Ids_file, 'rb' ) )
Summaries = cPickle.load( bz2.BZ2File( Summaries_file, 'rb' ) )

# <markdowncell>

# To make it easier to access the data, we convert here paper entries into [named tuples](http://docs.python.org/2/library/collections.html#collections.namedtuple). This will allow us to refer to fields by keyword, rather than index.

# <codecell>

from collections import namedtuple

paper = namedtuple( 'paper', ['title', 'authors', 'year', 'doi'] )

Summaries = {
    id : paper( *paper_info )
    for (id, paper_info) in Summaries.iteritems()
    }

# <codecell>

Summaries['23144668']

# <codecell>

Summaries['23144668'].title

# <headingcell level=2>

# Dataset statistics

# <markdowncell>

# Plotting relies on [matplotlib](http://matplotlib.org/), which you can download from [here](http://matplotlib.org/downloads.html) ([NumPy](http://www.numpy.org/) is also required, and can be downloaded [here](http://sourceforge.net/projects/numpy/files/NumPy/1.7.1/)).

# <codecell>

import matplotlib.pyplot as plt

# show plots inline within the notebook
%matplotlib inline
# set plots' resolution
plt.rcParams['savefig.dpi'] = 100

# <headingcell level=3>

# Papers per year

# <markdowncell>

# Here, we will get information on how many papers in the dataset were published per year.
# 
# We'll be using the [Counter](http://docs.python.org/2/library/collections.html#collections.Counter) class to determine the number of papers per year.

# <codecell>

paper_year = [ int(p.year) for p in Summaries.itervalues() ]

from collections import Counter

papers_per_year = sorted( Counter(paper_year).items() )
print 'Number of papers in the dataset per year for the past decade:'
print papers_per_year[-10:]

# <markdowncell>

# Filtering results, to obain only papers since 1950:

# <codecell>

papers_per_year  = [
    (y,count)
    for (y,count) in papers_per_year
    if y >= 1950
    ]

years     = [ y     for (y,count) in papers_per_year ]
nr_papers = [ count for (y,count) in papers_per_year ]

print 'Number of papers in the dataset published since 1950: %d.' % sum(nr_papers)

# <codecell>

plt.bar( left=years, height=nr_papers, width=1.0 )

plt.xlim(1950,2016)
plt.xlabel( 'year' )
plt.ylabel( 'number of papers' );

# <headingcell level=3>

# Papers per author

# <markdowncell>

# Here, we will obtain the distribution characterizing the number of papers published by an author.

# <codecell>

# flattening out of the list of lists of authors
authors_expanded = [
    auth
    for paper in Summaries.itervalues()
    for auth in paper.authors
    ]

nr_papers_by_author = Counter( authors_expanded )

# <codecell>

print 'There are %d authors in the dataset with distinct names.\n' % len(nr_papers_by_author)
print '15 authors with greatest number of papers:'
print sorted( nr_papers_by_author.items(), key=lambda i:i[1] )[-15:]

# <codecell>

plt.hist( x=nr_papers_by_author.values(), bins=range(51), histtype='step' )
plt.yscale('log')
plt.xlabel('number of papers authored')
plt.ylabel('number of authors');

# <headingcell level=3>

# Authors per paper

# <codecell>

plt.hist( x=[ len(p.authors) for p in Summaries.itervalues() ], bins=range(20), histtype='bar', align='left', normed=True )
plt.xlabel('number of authors in one paper')
plt.ylabel('fraction of papers')
plt.xlim(0,15);

# <headingcell level=3>

# Most frequently occurring words in paper titles

# <codecell>

# assemble list of words in paper titles, convert them to lowercase, and remove eventual trailing '.'
title_words = Counter([
    ( word if word[-1]!='.' else word[:-1] ).lower()
    for paper in Summaries.itervalues()
    for word in paper.title.split(' ')
    if word != ''
    ])

# <codecell>

print len(title_words), 'distinct words occur in the paper titles.\n'
print '30 most frequently occurring words:'
print sorted( title_words.items(), key=lambda i:i[1] )[-30:]

