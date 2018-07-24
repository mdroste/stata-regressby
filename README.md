
regressby
=================================

[Overview](#overview)
| [Installation](#installation)
| [Usage](#usage)
| [Benchmarks](#benchmarks)
| [To-Do](#todo)
| [Acknowledgements](#acknowledgements)
| [License](#license)

Flexible and hyper-fast grouped regressions in Stata

`version 0.5 24jul2018`


Overview
---------------------------------

regressby is a fast and efficient method to run grouped regressions; that is, it runs the same regression model on each subset partitioning your dataset. Functionally, this makes it very similar to the built-in -statsby- program, however, it runs between 10 and 100 times faster in most use cases. 


Installation
---------------------------------

There are two options for installing regressby.

1. The most recent version can be installed from Github with the following Stata commands:

```stata
local github "https://raw.githubusercontent.com"
net install regressby, from(`github'/mdroste/stata-regressby/master/)
```

2. A ZIP containing the program can be downloaded and manually placed on the user's adopath from Github.


Usage
---------------------------------

Concretely, consider an example where you have some variable y, x, and g. You want to regress y on x within each group indexed by a variable g.

The following two commands are equivalent:

``` 
regressby y x, by(byvars)
statsby, by(byvars) clear: reg y x	
```



Benchmarks
---------------------------------

Coming soon!
  
Todo
---------------------------------

The following items will be addressed soon:

- [ ] Finish off this readme.md
- [ ] Provide a help file
- [ ] Finish benchmarking and provide a script to validate results

A port of this program in C would yield a significant increase in performance; I have no plans to do that in the near future.


Acknowledgements
---------------------------------

This program is based off of internal code from the illustrious Michael Stepner's health inequality project. This program also benefited from contributions provided by the inimitable Dr. Wilbur Townsend, who helped elegantly generalize the code to allow for an arbitrary number of regressors. Finally, this program benefited greatly from the guidance and advice of Raj Chetty.


License
---------------------------------

regressby is [MIT-licensed](https://github.com/mdroste/stata-regressby/blob/master/LICENSE).

