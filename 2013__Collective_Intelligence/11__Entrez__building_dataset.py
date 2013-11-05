# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Building the _"evolution"_ research papers dataset

# <markdowncell>

# [Luís F. Simões](mailto:luis.simoes@vu.nl)<br>
# 2013-10-29 *(updated: 2013-11-05)*<div style="float: right">`Notebooks:` [next &rarr;](http://nbviewer.ipython.org/urls/raw.github.com/lfsimoes/VU/master/2013__Collective_Intelligence/12__inspecting_the_data.ipynb) &bull; [index &uarr;](https://github.com/lfsimoes/VU/tree/master/2013__Collective_Intelligence)</div><br><br>
# 
# *******

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
# </div><br><br>

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
# </div><br><br>

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

# <codecell>

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
# </div><br><br>

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

# ELink: Searching for related items in NCBI Entrez

# <markdowncell>

# <div style="float: right">
# `Entrez.elink()`: [Biopython Tutorial](http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc134), [module documentation](http://biopython.org/DIST/docs/api/Bio.Entrez-module.html#elink), [NCBI's E-utilities reference](http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ELink).
# </div><br><br>

# <markdowncell>

# To understand how to obtain paper citations with Entrez, we will first assemble a small set of PubMed IDs, and then query for their citations.
# To that end, we search here for papers published by [Chris Adami](http://adamilab.msu.edu/) in the [PLOS Computational Biology](http://www.ploscompbiol.org/) journal (as before, having also the word "evolution" in either the title or abstract):

# <codecell>

CA_search_term = search_term+'[TIAB] AND Adami C[AUTH] AND PLoS computational biology[JOUR]'
CA_ids = Entrez.read( Entrez.esearch( db="pubmed", term=CA_search_term ) )['IdList']
CA_ids

# <codecell>

CA_summ = {
    p['Id'] : ( p['Title'], p['AuthorList'], p['PubDate'][:4], p['FullJournalName'], p.get('DOI', '') )
    for p in Entrez.read( Entrez.esummary(db="pubmed", id=','.join( CA_ids )) )
    }
CA_summ

# <markdowncell>

# Because we restricted our search to papers in an open-access journal, you then follow their DOIs to freely access their PDFs at the journal's website: [10.1371/journal.pcbi.0040023](http://dx.doi.org/10.1371/journal.pcbi.0040023), [10.1371/journal.pcbi.1000948](http://dx.doi.org/10.1371/journal.pcbi.1000948), [10.1371/journal.pcbi.1002236](http://dx.doi.org/10.1371/journal.pcbi.1002236).

# <markdowncell>

# We will now issue calls to `Entrez.elink()` using these PubMed IDs, to retrieve the IDs of papers that cite them.
# The database from which the IDs will be retrieved is [PubMed Central](http://www.ncbi.nlm.nih.gov/pmc/), a free digital database of full-text scientific literature in the biomedical and life sciences.
# 
# You can, for instance, find [archived here](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2951343/), with the PubMed Central ID 2951343, the paper "Critical dynamics in the evolution of stochastic strategies for the iterated prisoner's dilemma", which as we saw above, has the PubMed ID 20949101.
# 
# A complete list of the kinds of links you can retrieve with `Entrez.elink()` can be found [here](http://eutils.ncbi.nlm.nih.gov/entrez/query/static/entrezlinks.html).

# <codecell>

CA_citing = {
    id : Entrez.read( Entrez.elink(
            cmd = "neighbor",               # ELink command mode: "neighbor", returns
                                            #     a set of UIDs in `db` linked to the input UIDs in `dbfrom`.
            dbfrom = "pubmed",              # Database containing the input UIDs: PubMed
            db = "pmc",                     # Database from which to retrieve UIDs: PubMed Central
            LinkName = "pubmed_pmc_refs",   # Name of the Entrez link to retrieve: "pubmed_pmc_refs", gets
                                            #     "Full-text articles in the PubMed Central Database that cite the current articles"
            from_uid = id                   # input UIDs
            ) )
    for id in CA_ids
    }

CA_citing['20949101']

# <markdowncell>

# We have in `CA_citing[paper_id][0]['LinkSetDb'][0]['Link']` the list of papers citing `paper_id`. To get it as just a list of ids, we can do

# <codecell>

cits = [ l['Id'] for l in CA_citing['20949101'][0]['LinkSetDb'][0]['Link'] ]
cits

# <markdowncell>

# However, one more step is needed, as what we have now are PubMed Central IDs, and not PubMed IDs. Their conversion can be achieved through an additional call to `Entrez.elink()`:

# <codecell>

cits_pm = Entrez.read( Entrez.elink( dbfrom="pmc", db="pubmed", LinkName="pmc_pubmed", from_uid=",".join(cits)) )
cits_pm

# <codecell>

ids_map = { pmc_id : link['Id'] for (pmc_id,link) in zip(cits_pm[0]['IdList'], cits_pm[0]['LinkSetDb'][0]['Link']) }
ids_map

# <markdowncell>

# And to check that indeed these are the PubMed IDs of papers citing the paper "Critical dynamics in the evolution of stochastic strategies for the iterated prisoner's dilemma", you can match PubMed Central's [list of citations for this paper](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2951343/citedby/) against the following results from a query to PubMed for their summaries:

# <codecell>

{   p['Id'] : ( p['Title'], p['AuthorList'], p['PubDate'][:4], p['FullJournalName'], p.get('DOI', '') )
    for p in Entrez.read( Entrez.esummary(db="pubmed", id=','.join( ids_map.values() )) )
    }

# <markdowncell>

# For comparison, you can match these with Google Scholar's [list of citations](http://scholar.google.com/scholar?cites=13285237943045850451) for the same paper.

# <headingcell level=3>

# Citations dataset

# <markdowncell>

# We have now seen all the steps required to assemble a dataset of citations to each of the papers in our dataset.

# <codecell>

Citations_file = search_term + '__Citations.pkl.bz2'
Citations = []

# <markdowncell>

# At least one server query will be issued per paper in `Ids`. Because NCBI allows for at most 3 queries per second (see [here](http://biopython.org/DIST/docs/api/Bio.Entrez-pysrc.html#_open)), this dataset will take a long time to assemble. Should you need to interrupt it for some reason, or the connection fail at some point, it is safe to just rerun the cell below until all data is collected.

# <codecell>

import httplib

if Citations == [] and os.path.exists( Citations_file ):
    Citations = cPickle.load( bz2.BZ2File( Citations_file, 'rb' ) )

if len(Citations) < len(Ids):
    
    i = len(Citations)
    checkpoint = len(Ids) / 10 + 1      # save to hard drive at every 10% of Ids fetched
    
    for pm_id in Ids[i:]:               # either starts from index 0, or resumes from where we previously left off
        
        while True:
            try:
                # query for papers archived in PubMed Central that cite the paper with PubMed ID `pm_id`
                c = Entrez.read( Entrez.elink( dbfrom = "pubmed", db="pmc", LinkName = "pubmed_pmc_refs", id=str(pm_id) ) )
                
                c = c[0]['LinkSetDb']
                if len(c) == 0:
                    # no citations found for the current paper
                    c = []
                else:
                    c = [ l['Id'] for l in c[0]['Link'] ]
                    
                    # convert citations from PubMed Central IDs to PubMed IDs
                    p = []
                    for start in xrange( 0, len(c), 200 ):
                        query_ids = ','.join( c[start : start+200] )
                        r = Entrez.read( Entrez.elink( dbfrom="pmc", db="pubmed", LinkName="pmc_pubmed", from_uid=query_ids ) )
                        # select the IDs. If no matching PubMed ID was found, [] is returned instead
                        p.extend( [] if r[0]['LinkSetDb']==[] else [ int(link['Id']) for link in r[0]['LinkSetDb'][0]['Link'] ] )
                    c = p
            
            except httplib.BadStatusLine:
                # Presumably, the server closed the connection before sending a valid response. Retry until we have the data.
                print 'r',
                continue
            break
        
        Citations.append( (pm_id, c) )
        i += 1
        print '\r at %.3f %%' % (100. * i/len(Ids)),
        
        if i % checkpoint == 0:
            print '\tsaving at checkpoint', i
            cPickle.dump( Citations, bz2.BZ2File( Citations_file, 'wb' ), protocol=pickle_protocol )
    
    print '\n done.'
    
    # Save Citations, as a dictionary indexed by Ids
    Citations = dict( Citations )
    
    cPickle.dump( Citations, bz2.BZ2File( Citations_file, 'wb' ), protocol=pickle_protocol )

# <markdowncell>

# To see that we have indeed obtained the data we expected, you can match the ids below (citations to the paper "Critical dynamics in the evolution of stochastic strategies for the iterated prisoner's dilemma"), with the ids listed at the end of last section.

# <codecell>

Citations[20949101]

# <headingcell level=2>

# Where do we go from here?

# <markdowncell>

# Running the code above generates multiple local files, containing the datasets we'll be working with. Loading them into memory is a matter of just issuing a call like<br>
# ``data = cPickle.load( bz2.BZ2File( data_file, 'rb' ) )``.
# 
# The Entrez module will therefore no longer be needed, unless you wish to extend your data processing with additional information retrieved from PubMed.
# 
# Should you be interested in looking at alternative ways to handle the data, have a look at the [sqlite3](http://docs.python.org/2/library/sqlite3.html) module included in Python's standard library, or [Pandas](http://pandas.pydata.org/), the Python Data Analysis Library.

