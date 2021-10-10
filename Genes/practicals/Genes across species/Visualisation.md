[Up to table of contents](README.md)

[Back to the previous page](How_much_of_the_genome_is_in_genes.md)

[Go to the next page](Scaling_up.md)

## A note on visualisation

It would be very nice to plot some of the values we've computed! For example - gene length versus
number of exons, or gene length versus number of transcripts, or genome size versus number of
genes, and so on.

Data visualisation is really beyond the scope of this tutorial. However, as a quick example, here
is how you might plot some of these values in R using [ggplot2](https://ggplot2.tidyverse.org):

```
library( RSQLite )
library( ggplot2 )
db = DBI::dbConnect( RSQLite::SQLite(), "genes.sqlite" )
data = dbGetQuery( db, "SELECT * FROM gene_summary_view" )
data$length = data$end - data$start + 1
p = (
   ggplot( data = data )
   + geom_point( aes( x = end - start + 1, y = average_number_of_exons, col = analysis ))
   + xlab( "Gene length" )
   + ylab( "Number of exons" )
)

ggsave( p, file = "gene_length_versus_number_of_exons.pdf", width = 8, height = 5 )
```

```
In python, you could instead try using  [matplotlib](https://matplotlib.org). I'd like to
include a matplotlib example here but don't have one - let me know if you write one!

[Go to the next page](Scaling_up.md)
