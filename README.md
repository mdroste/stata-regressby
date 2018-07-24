
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

regressby is a fast and efficient method to run grouped regressions; that is, it runs the same regression model on each subset partitioning your dataset and returns to you the coefficients. Functionally, this makes it very similar to the built-in -statsby- program, however, -regressby runs between 10 and 500 times faster than -statsby- in most use cases. These benefits are particularly large when there are many groups, when the number observations in each group is relatively small, and when the regression model only contains a few parameters.

Regressby supports robust and clustered standard errors, and also supports analytical weights.


Installation
---------------------------------

There are two options for installing regressby.

1. The most recent version can be installed from Github with the following Stata commands:

```stata
local github "https://raw.githubusercontent.com"
net install regressby, from(`github'/mdroste/stata-regressby/master/)
```

2. A ZIP containing the program can be downloaded and manually placed on the user's adopath from Github.


Motivating Example
---------------------------------

It is easiest to describe how regressby functions by way of example. Suppose you want to estimate a regression describing the relationship between a person's income, y, and their parent's income, x. Also suppose that you have a variable g that describes each person's place of birth (say, each value of g represents a county in the United States), and you would like to obtain slopes and intercepts separately for each county in the United States (as in Chetty and Hendren, 2014). 

You can accomplish this in one step by regressing y on a vector of dummy variables for each distinct value of g and a vector of interactions between these dummies and x. This is convenient, but potentially undesirable for a few reasons. For one, Stata isn't going to allow you to include more than 10,998 independent variables in your regression, so if you are interested in estimating group-specific slopes and intercepts, then you can only do this with fewer than 5,500 groups. Because there are only about 3,000 counties in the United States, it is feasible to perform this one-step estimation in Stata. However, it turns out that directly estimating thousands of parameters in an OLS regression simultaneously is quite slow in Stata.

An alternative estimatin strategy - and the only one that is feasible in a context where you have tens of thousands of distinct groups - is to estimate a separate regression of y on x for each distinct value of g. There are at least two easy ways to do this in Stata, either by manually iterating over values of groups or by using the built-in -statsby- function. However, both of these methods are also excruciatingly slow. Although it is very hard to beat Stata's built-in -regress- performance for a single regression, it can be very inefficient when running many identical regressions on subsets of your data.

The following two commands are equivalent:

``` 
regressby y x, by(byvars)
statsby, by(byvars) clear: reg y x	
```



Benchmarks
---------------------------------


![regressby benchmark](benchmarks/regressby_benchmark.png "regressby benchmark")
  
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

