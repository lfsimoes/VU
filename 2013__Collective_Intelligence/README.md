## Collective Intelligence ##
<a href="http://shop.oreilly.com/product/9780596529321.do">
<img align=right src="http://akamaicovers.oreilly.com/images/9780596529321/cat.gif"></a>

### About ###

This folder contains the code for generating the datasets used in the *Collective Intelligence*
course (FEW_X_401047_2013_120), along with examples of their usage.

The course uses the book 
"[Programming Collective Intelligence: Building Smart Web 2.0 Applications][CI_book]", by Toby Segaran (O'Reilly Media, 2007).

The provided datasets are meant to be used in implementations of the theory covered in the chapters on
"Making Recommendations" (Ch. 2), "Discovering Groups" (Ch. 3), "Searching and Ranking" (Ch. 4),
"Document Filtering" (Ch. 6) and "Finding Independent Features" (Ch. 10).


### IPython Notebooks ###
The code is provided in the form of IPython Notebooks, each of which is also available as a .py file export.

The links below visualize the the IPython Notebooks contained in this folder through the [IPython Notebook Viewer](http://nbviewer.ipython.org/).

1. [Building the "evolution" research papers dataset][Notebook_11]
2. [Inspecting the dataset][Notebook_12]
3. [Network analysis][Notebook_13]
4. [Text mining][Notebook_14]

An extra dataset, containing papers' keywords, was also provided, as a form of ground truth of paper categorization.
Its construction used the code in the [Extra__Keywords][Notebook_KW] notebook.


### Follow-ups ###

The VU Amsterdam's [Information Retrieval][VUIR] course has since 2014 provided exercises derived from these notebooks.

See the repositories for:

* [2014/2015][IR1415] (by [Paul Groth][pgroth], with assignments by J.E. Hoeksema)
* [2015/2016][IR1516], [2016/2017][IR1617], [fall 2017][IR1718], [fall 2018][IRf18] (by [Tobias Kuhn][tkuhn]) &ndash; Converted to the latest Jupyter version and Python 3



[Notebook_11]: http://nbviewer.ipython.org/urls/raw.github.com/lfsimoes/VU/master/2013__Collective_Intelligence/11__Entrez__building_dataset.ipynb
[Notebook_12]: http://nbviewer.ipython.org/urls/raw.github.com/lfsimoes/VU/master/2013__Collective_Intelligence/12__inspecting_the_data.ipynb
[Notebook_13]: http://nbviewer.ipython.org/urls/raw.github.com/lfsimoes/VU/master/2013__Collective_Intelligence/13__network_analysis.ipynb
[Notebook_14]: http://nbviewer.ipython.org/urls/raw.github.com/lfsimoes/VU/master/2013__Collective_Intelligence/14__text_mining.ipynb
[Notebook_KW]: http://nbviewer.ipython.org/urls/raw.github.com/lfsimoes/VU/master/2013__Collective_Intelligence/Extra__Keywords.ipynb

[CI_book]: http://shop.oreilly.com/product/9780596529321.do

[VUIR]: https://github.com/VUInformationRetrieval
[IR1415]: https://github.com/VUInformationRetrieval/IR2014_2015
[IR1516]: https://github.com/VUInformationRetrieval/IR2015_2016
[IR1617]: https://github.com/VUInformationRetrieval/IR2016_2017
[IR1718]: https://github.com/VUInformationRetrieval/IR2017_2018
[IRf18]: https://github.com/VUInformationRetrieval/IR2018
[pgroth]: https://github.com/pgroth
[tkuhn]: https://github.com/tkuhn
