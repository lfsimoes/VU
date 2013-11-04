# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Building the _"evolution"_ research papers dataset

# <markdowncell>

# [Luís F. Simões](mailto:luis.simoes@vu.nl)<br>
# 2013-10-29 (*updated: 2013-11-02*)<br><br><br>

# <markdowncell>

# The [Entrez](http://biopython.org/DIST/docs/api/Bio.Entrez-module.html) module, a part of the [Biopython](http://biopython.org/) library, will be used to interface with [PubMed](http://www.ncbi.nlm.nih.gov/pubmed).<br>
# You can download Biopython from [here](http://biopython.org/wiki/Download).
# 
# In this notebook we will be covering several of the steps taken in the [Biopython Tutorial](http://biopython.org/DIST/docs/tutorial/Tutorial.html), specifically in [Chapter 9  Accessing NCBI’s Entrez databases](http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc109).

# <codecell>

from Bio import Entrez

# NCBI requires you to set your email address to make use of NCBI's E-utilities
Entrez.email = "Your.Name.Here@example.org"

# <markdowncell>

# The datasets will be saved as serialized Python objects, compressed with bzip2.
# Saving/loading them will therefore require the [cPickle](http://docs.python.org/2/library/pickle.html) and [bz2](http://docs.python.org/2/library/bz2.html) modules.

# <codecell>

import cPickle, bz2
import os

# data format used in Python object serialization
pickle_protocol = 2

# <headingcell level=2>

# EInfo: Obtaining information about the Entrez databases

# <markdowncell>

# <div style="float: right">
# `Entrez.einfo()`: [Biopython Tutorial](http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc111), [module documentation](http://biopython.org/DIST/docs/api/Bio.Entrez-module.html#einfo), [NCBI's E-utilities reference](http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EInfo).
# </div>

# <codecell>

# accessing extended information about the PubMed database
pubmed = Entrez.read( Entrez.einfo(db="pubmed"), validate=False )[u'DbInfo']

# list of possible search fields for use with ESearch:
search_fields = { f['Name']:f['Description'] for f in pubmed["FieldList"] }

# <markdowncell>

# In search_fields, we find 'TIAB' ('Free text associated with Abstract/Title') as a possible search field to use in searches.

# <codecell>

search_fields

# <headingcell level=2>

# ESearch: Searching the Entrez databases

# <markdowncell>

# <div style="float: right">
# `Entrez.esearch()`: [Biopython Tutorial](http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc112), [module documentation](http://biopython.org/DIST/docs/api/Bio.Entrez-module.html#esearch), [NCBI's E-utilities reference](http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch).
# </div>

# <markdowncell>

# To have a look at the kind of data we get when searching the database, we'll perform a search for papers authored by Haasdijk:

# <codecell>

example_authors = ['Haasdijk E']
example_search = Entrez.read( Entrez.esearch( db="pubmed", term=' AND '.join([a+'[AUTH]' for a in example_authors]) ) )
example_search

# <markdowncell>

# Note how the result being produced is not in Python's native string format:

# <codecell>

type( example_search['IdList'][0] )

# <markdowncell>

# The part of the query's result we are most interested in is accessible through

# <codecell>

example_ids = [ int(id) for id in example_search['IdList'] ]
print example_ids

# <headingcell level=3>

# PubMed IDs dataset

# <markdowncell>

# We will now assemble a dataset comprised of research articles containing the keyword "evolution", in either their titles or abstracts.

# <codecell>

search_term = 'evolution'

Ids_file = search_term + '__Ids.pkl.bz2'

# <codecell>

if os.path.exists( Ids_file ):
    Ids = cPickle.load( bz2.BZ2File( Ids_file, 'rb' ) )
else:
    # determine the number of hits for the search term
    search = Entrez.read( Entrez.esearch( db="pubmed", term=search_term+'[TIAB]', retmax=0 ) )
    total = int( search['Count'] )
    
    # `Ids` will be incrementally assembled, by performing multiple queries,
    # each returning at most `retrieve_per_query` entries.
    Ids_str = []
    retrieve_per_query = 10000
    
    for start in xrange( 0, total, retrieve_per_query ):
        print 'Fetching IDs of results [%d,%d]' % ( start, start+retrieve_per_query )
        s = Entrez.read( Entrez.esearch( db="pubmed", term=search_term+'[TIAB]', retstart=start, retmax=retrieve_per_query ) )
        Ids_str.extend( s[ u'IdList' ] )
    
    # convert Ids to integers (and ensure that the conversion is reversible)
    Ids = [ int(id) for id in Ids_str ]
    
    for (id_str, id_int) in zip(Ids_str, Ids):
        if str(id_int) != id_str:
            raise Exception('Conversion of PubMed ID %s from string to integer it not reversible.' % id_str )
    
    # Save list of Ids
    cPickle.dump( Ids, bz2.BZ2File( Ids_file, 'wb' ), protocol=pickle_protocol )
    
total = len( Ids )
print '%d documents contain the search term "%s".' % ( total, search_term )

# <markdowncell>

# Taking a look at what we just retrieved, here are the last 5 elements of the `Ids` list:

# <codecell>

Ids[:5]

# <headingcell level=2>

# ESummary: Retrieving summaries from primary IDs

# <markdowncell>

# <div style="float: right">
# `Entrez.esummary()`: [Biopython Tutorial](http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc114), [module documentation](http://biopython.org/DIST/docs/api/Bio.Entrez-module.html#esummary), [NCBI's E-utilities reference](http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESummary).
# </div>

# <markdowncell>

# To have a look at the kind of metadata we get from a call to `Entrez.esummary()`, we now fetch the summary of one of Haasdijk's papers (using one of the PubMed IDs we obtained in the [previous section](#ESearch:-Searching-the-Entrez-databases)):

# <codecell>

example_paper = Entrez.read( Entrez.esummary(db="pubmed", id='23144668') )[0]

def print_summary( p ):
    for k,v in p.items():
        print k
        print '\t',v

print_summary(example_paper)

# <markdowncell>

# For now, we'll keep just some basic information for each paper: title, list of authors, publication year, and [DOI](https://en.wikipedia.org/wiki/Digital_object_identifier).
# 
# In case you are not familiar with the DOI system, know that the paper above can be accessed through the link [http://dx.doi.org/10.1007/s12065-012-0071-x](http://dx.doi.org/10.1007/s12065-012-0071-x) (which is `http://dx.doi.org/` followed by the paper's DOI).

# <codecell>

( example_paper['Title'], example_paper['AuthorList'], int(example_paper['PubDate'][:4]), example_paper['DOI'] )

# <headingcell level=3>

# Summaries dataset

# <markdowncell>

# We are now ready to assemble a dataset containing the summaries of all the paper `Ids` we previously fetched.
# 
# To reduce the memory footprint, and to ensure the saved datasets won't depend on Biopython being installed to be properly loaded, values returned by `Entrez.read()` will be converted to their corresponding native Python types. We start by defining a function for helping with the conversion of strings:

# <codecell>

def conv_str( s ):
    """
    The Entrez parser returns strings as either instances of
    `Bio.Entrez.Parser.StringElement` (a subclass of `str`),
    or `Bio.Entrez.Parser.UnicodeElement` (a subclass of `unicode`).
    Such strings are here converted to Python's native formats.
    
    See: http://biopython.org/DIST/docs/api/Bio.Entrez.Parser-pysrc.html
    """
    return unicode(s) if isinstance(s, unicode) else str(s)

# <codecell>

Summaries_file = search_term + '__Summaries.pkl.bz2'

# <codecell>

if os.path.exists( Summaries_file ):
    Summaries = cPickle.load( bz2.BZ2File( Summaries_file, 'rb' ) )
else:
    # `Summaries` will be incrementally assembled, by performing multiple queries,
    # each returning at most `retrieve_per_query` entries.
    Summaries = []
    retrieve_per_query = 200
    
    print 'Fetching Summaries of results: ',
    for start in xrange( 0, len(Ids), retrieve_per_query ):
        print start,
        
        # build comma separated string with the ids at indexes [start, start+retrieve_per_query)
        query_ids = ','.join( [ str(id) for id in Ids[ start : start+retrieve_per_query ] ] )
        
        s = Entrez.read( Entrez.esummary( db="pubmed", id=query_ids ) )
        
        # out of the retrieved data, we will keep only a tuple (title, authors, year, DOI), associated with the paper's id.
        # (all values converted to native Python formats)
        f = [
            ( int( p['Id'] ), (
                conv_str( p['Title'] ),
                [ conv_str(a) for a in p['AuthorList'] ],
                int( p['PubDate'][:4] ),                # keeps just the publication year
                conv_str( p.get('DOI', '') )            # papers for which no DOI is available get an empty string in their place
                ) )
            for p in s
            ]
        Summaries.extend( f )
    
    # Save Summaries, as a dictionary indexed by Ids
    Summaries = dict( Summaries )
    
    cPickle.dump( Summaries, bz2.BZ2File( Summaries_file, 'wb' ), protocol=pickle_protocol )

# <markdowncell>

# Let us take a look at the first 3 retrieved summaries:

# <codecell>

{ id : Summaries[id] for id in Ids[:3] }

# <headingcell level=2>

# Where do we go from here?

# <markdowncell>

# Running the code above generates multiple local files, containing the datasets we'll be working with. Loading them into memory is a matter of just issuing a call like<br>
# ``data = cPickle.load( bz2.BZ2File( data_file, 'rb' ) )``.
# 
# The Entrez module will therefore no longer be needed, unless you wish to extend your data processing with additional information retrieved from PubMed.
# 
# Should you be interested in looking at alternative ways to handle the data, have a look at the [sqlite3](http://docs.python.org/2/library/sqlite3.html) module included in Python's standard library, or [Pandas](http://pandas.pydata.org/), the Python Data Analysis Library.

