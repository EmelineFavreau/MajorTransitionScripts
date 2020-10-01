Learning about unit testing

Why doing it?
Because mistakes happen all the times, like in 1999 when NASA sent a $125 million satellite to burn off because one of the software they used was set to newtons and not pounds of force.
So to save ourselves some time, and to write code that is easier to use (and less awkward to feel about), testing the code in units is the way forward.

How to start?
- write any test, even a small one like ifelse statement
- run the tests often - can you automatize it? can you make it part of the script template?
- it's better to run a short test (on a chunk of code) than a long test (on the whole software/pipeline)
- it's better to have a good coverage of the script
- document the test with comments, toydataset that are realistic to your real data (e.g. with as many null values as the real data), etc
- use tools: pytest (recommended) and unittest comes with standard Python library; R unit-tests: testthat, tinytest, svUnit

Types of testing:
- smoke testing: initial checks for anything badly broken. 
- unit test: targeted to one code chunk, like using toydataset on one R function to check output class
- integration testing: testing communication between code chunks. Very important when writing a pipeline with different people (otherwise you burn the expensive satellite, see above)
- system tests: from end-to-end programs, they take a long time
- regression testing: compare the output of the code after it has been changed. Very important in team working.
- runtime testing: they are part of the program and alert if something is off, like negative population level

Special cases:
- stochastic code (e.g. Monte Carlo methods) will likely fail during unit testing, because of different outputs. Using ransom number seeds or exploring distribution of results can help.
- sometimes the data need eyeballing (human curation).
- floating point issues in Python

More reading:
- https://the-turing-way.netlify.app/reproducible-research/testing.html
- R testthat: https://testthat.r-lib.org/
- testing in R: https://r-pkgs.org/tests.html