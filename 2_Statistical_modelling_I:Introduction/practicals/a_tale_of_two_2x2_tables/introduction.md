# A tale of two 2x2 tables

Welcome to (possibly) the world's first *choose your own statistical adventure*!

In this practical we are going to study the following two two-by-two tables, which come from two
published papers about genetic variants that affect malaria susceptibility:
<https://doi.org/10.1038/s41467-019-13480-z> and
<https://science.sciencemag.org/content/348/6235/711>.

```R
                       non-O         O                                                   rs60822373     rs60822373
                       blood group   blood group                                         G allele       C allele
            controls:  3420          3233                        unexposed populations:  1965           1
severe malaria cases:  3925          2738                  malaria-exposed populations:  707            17
```

The left table purports to show evidence that O blood group is associated with protection against
severe malaria. The right table purports to show evidence that rs60822373 has evolved under
selection from malaria (that has kept it at higher frequency in malaria-exposed populations).

If you run a statistical test on either table you will see that there are indeed highly
statistically significant differences in frequencies between the two rows of each table, so these
differences are not just due to sampling.

But *something is wrong with one of these tables*.  Your mission is to find out what!

*QN*.  Maybe you can already see what is wrong?

To get started let's load the tables into R:
```R
tables = list(
    casecontrol = matrix(
        c( 3420, 3233, 3925, 2738 ),
        byrow = T, ncol = 2,
        dimnames = list(
            c( "non-O", "O" ),
            c( "controls", "cases" )
        )
    ),
    exposure = matrix(
        c( 1965, 1, 707, 17 ),
        byrow = T, ncol = 2,
        dimnames = list(
            c( "G", "C" ),
            c( "unexposed", "exposed" )
        )
    )
)
```

[I'm suspicious about the data.](suspicious_data.md)
[I'm suspicious about the statistical test.](suspicious_test.md)
