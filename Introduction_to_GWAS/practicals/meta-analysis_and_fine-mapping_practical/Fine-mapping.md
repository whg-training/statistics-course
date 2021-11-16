[Up to the table of contents](README.md) - [Back to the meta-analysis sectioin](meta-analysis.md))

### Using `FINEMAP` to fine-map associations

If you followed the [meta-analysis section](Meta-analysis.md) you should now have a dataframe with meta-analysis results.  It looks something like this:

    > meta_analysis
    # A tibble: 730 x 16
       rsid             chromosome position allele1 allele2 study1.beta study1.se   study1.P log10_study1.BF study2.beta study2.se study2.P meta.beta meta.se   meta.P log10_meta.BF
       <chr>                 <dbl>    <dbl> <chr>   <chr>         <dbl>     <dbl>      <dbl>           <dbl>       <dbl>     <dbl>    <dbl>     <dbl>   <dbl>    <dbl>         <dbl>
     1 19:49149607:G:A          19 49149607 G       A            0.202     0.228  0.376              -0.0925     -0.146      0.295  0.310      0.0176  0.0894 8.44e- 1       -0.0497
     2 19:49149711:G:A          19 49149711 G       A           -0.0183    0.202  0.928               0.0195      0.239      0.261  0.180      0.0245  0.0894 7.84e- 1       -0.0668
     3 19:49149937:T:G          19 49149937 T       G           -0.357     0.0785 0.00000540          4.90       -0.264      0.101  0.00463   -0.669   0.0894 7.45e-14       12.3   
     4 19:49150048:G:T          19 49150048 G       T           -0.332     0.0788 0.0000257           4.26       -0.271      0.102  0.00387   -0.636   0.0894 1.12e-12       11.2   
     5 19:49150130:G:T          19 49150130 G       T           -0.196     0.110  0.0756              0.947      -0.204      0.143  0.0767    -0.209   0.0894 1.96e- 2        1.52  
     6 19:49150233:G:A          19 49150233 G       A            0.133     0.134  0.321              -0.137      -0.0794     0.173  0.323      0.0381  0.0894 6.70e- 1       -0.0958
     7 19:49151312:C:T          19 49151312 C       T            0.238     0.145  0.101              -0.137       0.101      0.187  0.295      0.114   0.0894 2.03e- 1       -0.176 
     8 19:49151323:T:G          19 49151323 T       G            0.202     0.228  0.376              -0.0925     -0.146      0.295  0.310      0.0176  0.0894 8.44e- 1       -0.0497
     9 19:49151460:C:CA         19 49151460 C       CA          -0.283     0.153  0.0636              0.936      -0.218      0.197  0.134     -0.142   0.0894 1.12e- 1        0.830 
    10 19:49152118:C:T          19 49152118 C       T           -0.609     0.328  0.0632              0.634      -0.604      0.423  0.0767    -0.0723  0.0894 4.19e- 1        0.319 
