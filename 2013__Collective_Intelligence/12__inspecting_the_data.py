# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Inspecting the dataset

# <markdowncell>

# [Luís F. Simões](mailto:luis.simoes@vu.nl), 2013-11-01<br><br><br>

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

#Ids = cPickle.load( bz2.BZ2File( Ids_file, 'rb' ) )
Summaries = cPickle.load( bz2.BZ2File( Summaries_file, 'rb' ) )

# <markdowncell>

# To make it easier to access the data, we convert here paper entries into [named tuples](http://docs.python.org/2/library/collections.html#collections.namedtuple). This will allow us to refer to fields by keyword, rather than index.

# <codecell>

from collections import namedtuple

paper = namedtuple( 'paper', ['title', 'authors', 'year', 'doi'] )

for (id, paper_info) in Summaries.iteritems():
    Summaries[id] = paper( *paper_info )

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

# <markdowncell>

# Creating a bar plot to visualize the results (using [matplotlib.pyplot.bar](http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.bar)):

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

# <markdowncell>

# Creating a histogram to visualize the results (using [matplotlib.pyplot.hist](http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.hist)):

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
    if word != ''     # the split by spaces generates empty strings when consecutive spaces occur in the title; this discards them
    ])

# <codecell>

print len(title_words), 'distinct words occur in the paper titles.\n'
print '50 most frequently occurring words:'
print sorted( title_words.items(), key=lambda i:i[1] )[-50:]

# <headingcell level=3>

# Tracing the evolution of research topics

# <markdowncell>

# We will be measuring here how often did certain keywords occur in paper titles across different years.
# Remember though, that we are dealing with a biased dataset: we will be measuring frequencies of keywords among titles that were already preselected (for occurrence of "evolution" in their titles or abstracts).

# <markdowncell>

# Grouping paper titles by publication year, while discarding papers published before 1950 (using [collections.defaultdict](http://docs.python.org/2/library/collections.html#collections.defaultdict)):

# <codecell>

from collections import defaultdict

# `papers_per_year` will be a dictionary indexing year to titles of papers published in that year.
# As a defaultdict (of lists), actions on non-existing keys cause them to be added to the dictionary, with an empty list as value.
papers_per_year = defaultdict(list)

for p in Summaries.itervalues():
    if int(p.year) >= 1950:
        papers_per_year[ int(p.year) ].append( p.title )

years = sorted( papers_per_year.keys() )

# <markdowncell>

# Checking the outcome: first 5 paper titles in 1950

# <codecell>

papers_per_year[1950][:5]

# <markdowncell>

# Measuring fraction of papers in each year having the different keywords in their titles.<br>
# Beware that the way we are doing things here, a keyword like 'gene' matches also things like 'genetics', or even 'genealogy'.

# <codecell>

keywords = ['DNA', 'gene', 'genom', 'simulation', 'model']
#keywords = ['disease', 'virus', 'cancer', 'HIV']
#keywords = ['adaptation', 'speciation', 'extinction', 'mutation', 'variation']

keywords_per_year = {}

for kw in keywords:
    keywords_per_year[kw] = [
        sum([ 1 for title in papers_per_year[y] if kw in title ]) / float(len(papers_per_year[y]))
        for y in years
        ]

# <markdowncell>

# Listing a keyword's usage over a period of time:

# <codecell>

keywords_per_year['gene'][ years.index(2005) : years.index(2013)+1 ]

# <markdowncell>

# Creating a plot to visualize the results (using [matplotlib.pyplot.plot](http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.plot)):

# <codecell>

for kw in keywords:
    plt.plot( years, keywords_per_year[kw], label=kw )

plt.legend( loc='upper left', frameon=False )
plt.xlabel('year')
plt.ylabel('fraction of papers with keyword in title')
plt.xlim(1950,2013);

# <headingcell level=3>

# *Authored*, and *Coauthors* mappings

# <markdowncell>

# `Summaries` maps paper ids to paper summaries. Let us start by building a mapping from authors, to ids of papers they authored:

# <codecell>

papers_of_author = defaultdict(list)

for id,p in Summaries.iteritems():
    for a in p.authors:
        papers_of_author[ a ].append( id )

# <codecell>

papers_of_author['Eiben AE']

# <codecell>

{ pid : Summaries[pid] for pid in papers_of_author['Eiben AE'] }

# <markdowncell>

# From the `papers_of_author` mapping, we will now build a mapping from authors, to the set of co-authors they have published with:

# <codecell>

coauthors = defaultdict(set)

for a,pids in papers_of_author.iteritems():
    for id in pids:
        for ca in Summaries[id].authors:
            if ca != a:
                coauthors[ a ].add( ca )

# <codecell>

coauthors['Eiben AE']

# <markdowncell>

# With this data in hand, we can now plot the distribution showing the number of collaborators a scientist has published with in its full publication record:

# <codecell>

plt.hist( x=[ len(ca) for ca in coauthors.itervalues() ], bins=range(52), histtype='bar', align='left', normed=True )
plt.xlabel('number of collaborators')
plt.ylabel('fraction of scientists')
plt.xlim(0,50);

