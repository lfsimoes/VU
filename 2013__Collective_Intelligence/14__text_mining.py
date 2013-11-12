# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Text mining

# <markdowncell>

# [Luís F. Simões](mailto:luis.simoes@vu.nl)<br>
# 2013-11-08 *(updated: 2013-11-12)*<div style="float: right">`Notebooks:` [previous &larr;](http://nbviewer.ipython.org/urls/raw.github.com/lfsimoes/VU/master/2013__Collective_Intelligence/13__network_analysis.ipynb) &bull; [index &uarr;](https://github.com/lfsimoes/VU/tree/master/2013__Collective_Intelligence)</div><br><br>
# 
# *******

# <markdowncell>

# This notebook provides some basic examples of how to process textual data.
# 
# You have in the [Natural Language Toolkit (NLTK)](http://nltk.org/), and [CLiPS Pattern](http://www.clips.ua.ac.be/pages/pattern), sophisticated libraries for natural language processing with Python, which you are free to deploy over this data.

# <headingcell level=2>

# Loading modules & datasets

# <codecell>

import cPickle, bz2
from collections import *

import numpy as np

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

# [Bag-of-words model](https://en.wikipedia.org/wiki/Bag-of-words_model)

# <markdowncell>

# Here we'll find out how to represent text by vectors of word frequencies.
# 
# Given a string (or several), the function below will remove a few eventual characters from it, convert all characters to lower case, split it by spaces, and then return a counter showing the number of times each word occurs.
# 
# Should you prefer to bring heavier weaponry into this same task, have a look at Python's [regular expressions](http://docs.python.org/2/library/re.html) module.
# 
# > *Implementation note*: [list comprehension](http://docs.python.org/2/tutorial/datastructures.html#list-comprehensions) (additional examples [here](http://www.pythonforbeginners.com/lists/list-comprehensions-in-python/)) has been used extensively in these examples. In the function below, however, you'll note that the same syntax is used, but the square brackets are replaced by parentheses. Those are [generator expressions](http://www.python.org/dev/peps/pep-0289/), which create objects meant for a single iteration over the resulting sequence. Given the volume of data we are handling, such expressions can provide very significant savings in memory consumption, and also computation time.

# <codecell>

def word_freq( text ):
    """
    Splits text into words, and counts their occurrence frequencies.
    """
    if isinstance(text, basestring):
        # if a string (or unicode string) was sent, put it into a single element tuple
        text = (text,)
    
    # performs character substitutions/removals, converts to lower case, and splits by spaces
    text_clean = (
        txt.replace('/',' ').replace(',',' ').replace('.',' ').replace('(','').replace(')','').lower().split(' ')
        for txt in text
        )
    
    return Counter((
        word
        for txt  in text_clean
        for word in txt
        if word != ''   # the split by spaces generates empty strings when consecutive spaces occur in the title; this discards them
        ))

# <markdowncell>

# This function can accept a single string, but it can also calculate its frequencies over a sequence of strings (no need to concatenate them first).

# <codecell>

word_freq( 'this and that and also something else' )

# <codecell>

word_freq( ['here is a string', 'and here is another'] )

# <markdowncell>

# We now apply it to two papers' title and abstract strings.

# <codecell>

[   ( Summaries[7466396].title, Abstracts[7466396] ),
    ( Summaries[8316296].title, Abstracts[8316296] )
    ]

# <codecell>

p1 = word_freq( [Summaries[7466396].title, Abstracts[7466396]] )
p1

# <codecell>

p2 = word_freq( Summaries[8316296].title + ' ' + Abstracts[8316296] )
p2

# <markdowncell>

# Now, we use set operations to intersect the sets of words occurring in each paper, and then sort the results by decreasing order of total usage.

# <codecell>

shared = sorted( set( p1.keys() ) & set( p2.keys() ), key=lambda i:p1[i]+p2[i], reverse=True )
print shared

# <markdowncell>

# We'll use these words to create our vectors summarizing each paper's contents.

# <codecell>

word_vects = [
    [ freqs.get(word,0) for word in shared ]
    for freqs in [p1,p2]
    ]
word_vects

# <markdowncell>

# In this representation, we know that the word 'and' (the third one in the `shared` list) occurs 3 times in the first document, and 6 in the second.

# <markdowncell>

# Using [Numpy](http://www.numpy.org/) ([introductory reference](http://docs.scipy.org/doc/numpy/user/basics.html), [tutorial](http://wiki.scipy.org/Tentative_NumPy_Tutorial)) you can convert this from a list of lists into a matrix form:

# <codecell>

word_vects = np.array( word_vects )
word_vects

# <markdowncell>

# Using Numpy will allow you to perform many operations over whole matrices/vectors, with a compact syntax and very high performance (Numpy relies on C libraries for many of its operations).
# 
# Below we normalize each paper's vector:

# <codecell>

rows_sum = word_vects.sum( axis=1 )
rows_sum

# <codecell>

rows_sum.reshape(-1,1).astype( np.float )

# <codecell>

n = word_vects / rows_sum.reshape(-1,1).astype( np.float )
n

# <markdowncell>

# We can obtain just one column as such:

# <codecell>

n[ : , shared.index('strategy') ]

# <headingcell level=3>

# Simple word clouds

# <markdowncell>

# In this section, we will be printing word clouds representing words' occurrence frequencies in a corpus of text. We will use IPython's Rich Display System ([examples here](http://nbviewer.ipython.org/url/github.com/ipython/ipython/raw/master/examples/notebooks/Part%205%20-%20Rich%20Display%20System.ipynb)) to have code cells producing HTML output. This section also provides a more elaborate example of how to use Numpy.

# <codecell>

from IPython.display import display, HTML

def word_cloud( words, frequency ):
    """
    Given a sequence of `words`, and a matching sequence of `frequency` values (scaled in [0,1]),
    displays a word cloud with word sizes scaled proportionally to their frequencies.
    """
    max_font_size = 40.
    min_font_size = 8.

    wc = HTML( ' '.join(
        '<span style="font-size:%.1fpx;">&nbsp;%s&nbsp;</span>' % ( min_font_size + f*(max_font_size-min_font_size), w )
        for w,f in zip(words, frequency)
        ) )
    display(wc)

# <markdowncell>

# The corpus of text we will be using as an example is made up of the titles and abstracts of papers containing the word "extinct" in their abstracts (words such as "extinction", "Extinction", "[pseudoextinction](https://en.wikipedia.org/wiki/Pseudoextinction)" and "[coextinction](https://en.wikipedia.org/wiki/Coextinction)" will also be caught up in our net).
# 
# Note how a generator expression is used to create an iterator over all strings matching our criterion (which is then sent as argument to `word_freq()`). This way, we never actually build in memory the whole list of strings in our corpus. Of course, this only works here because `word_freq()` was coded in such a way that it can accept iterators as arguments.

# <codecell>

extinction_corpus_freqs = word_freq(
    element
    for (id,abstr) in Abstracts.iteritems()
    for element in (Summaries[id].title, abstr)
    if 'extinct' in abstr or 'Extinct' in abstr
    )

# <markdowncell>

# We now sort by decreasing order of word frequency, and show the 10 most frequently occurring words.

# <codecell>

extinction = sorted( extinction_corpus_freqs.items(), key=lambda i:i[1], reverse=True )
extinction[:10]

# <markdowncell>

# Using Python's [zip](http://docs.python.org/2/library/functions.html#zip) built-in function (with an [argument list unpacked](http://docs.python.org/2/tutorial/controlflow.html#unpacking-argument-lists) through `*`) we can easily decompose a list of tuples into multiple lists containing each just the elements at a certain position.

# <codecell>

extinction_w, extinction_freq = zip( *extinction )

# <markdowncell>

# Here are then the 20 most frequently occurring words, and their frequencies:

# <codecell>

print extinction_w[:20]

# <codecell>

print extinction_freq[:20]

# <markdowncell>

# We can now easily convert the list of word frequencies into a numpy array, for easier handling.

# <codecell>

extinction_freq = np.array( extinction_freq )

print 'Array shape:', extinction_freq.shape
print 'First 100 word frequencies:'
extinction_freq[:100]

# <markdowncell>

# If we were to scale word sizes based on these frequencies, we'd see the first few words, and little else. Because the most frequently occurring words are used so much more often than the rest, we can't really scale based on them. We'll therefore impose a frequency cutoff: words with frequencies greater than the cutoff will be brought down to it. We will set here the cutoff to the 99.75th percentile of the data (with [numpy.percentile](http://docs.scipy.org/doc/numpy/reference/generated/numpy.percentile.html)).

# <codecell>

freq_cutoff = np.percentile( extinction_freq, 99.75 )
freq_cutoff

# <codecell>

extinction_freq = np.clip( extinction_freq, a_min=0, a_max=freq_cutoff )
extinction_freq[:100]

# <markdowncell>

# We don't want to print out over 20 thousand words. So, we'll impose a maximum number of printed words, and proceed with scaling each word's frequency in [0,1] according to the range of values observed in the set of words we'll be printing.

# <codecell>

max_words = 250

# <codecell>

extinction_selection = extinction_freq[:max_words]

# <codecell>

m, M = extinction_selection.min(), extinction_selection.max()

extinction_selection = ( extinction_selection - m ) / (M - m)
extinction_selection[:100]

# <markdowncell>

# And, having concluded all the required intermediate steps, we can finally print our word cloud:

# <codecell>

word_cloud( extinction_w[:max_words], extinction_selection )

# <markdowncell>

# And there you have it. The collective intelligence of papers about extinction, tells us that:
# 
# > The of and in to a that is for evolution with are we from by as on species this extinction an have be or these which evolutionary at between population can has not their was may more been were model extinct but it populations.

# <markdowncell>

# Now that you understand all the calculations we perform with words' frequencies, we can condense `word_cloud()` so it will carry out all these operations, and enable it to take as input a Counter object returned from a call to `word_freq()`.

# <codecell>

from IPython.display import display, HTML

def word_cloud( word_frequencies, max_words=250, cutoff_percentile=90, font_size_range=(8,40) ):
    """
    Given a Counter object (such as the one returned by `word_freq()`), display a word cloud
    with the `max_words` most frequently occurring words.
    """
    words, freq = zip( *sorted( word_frequencies.items(), key=lambda i:i[1], reverse=True ) )
    
    freq = np.array( freq[:max_words] )
    freq = np.clip( freq, a_min=0, a_max=np.percentile(freq, cutoff_percentile) )
    
    # because freq was sorted by descending order of word frequencies, we can obtain its range of values at the extremes
    m = freq[-1]
    M = freq[0]
    freq = ( freq - m ) / (M - m)
    
    min_fs, max_fs = font_size_range
    
    wc = HTML( ' '.join(
        '<span style="font-size:%.1fpx;">&nbsp;%s&nbsp;</span>' % ( min_fs + f*(max_fs - min_fs), w )
        for w,f in zip(words[:max_words], freq[:max_words])
        ) )
    display(wc)

# <markdowncell>

# And as an example, here is the word cloud generated from all abstracts having the word "simulat" (simulate, simulation, ...) in them.

# <codecell>

word_cloud( word_freq( a for a in Abstracts.itervalues() if 'simulat' in a ) )

# <markdowncell>

# Many words in these word clouds are obviously not good descriptors of the papers they appear in (the, of, and...). This is where you now come in with [tf–idf](https://en.wikipedia.org/wiki/Tf%E2%80%93idf)...

