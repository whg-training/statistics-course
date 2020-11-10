n = 1
data = read.table( sprintf( "diagram%d_simulated_data.tsv", n ), hea=T, as.is=T, sep ="\t" )


# Fit model without B
l = lm( C ~ A, data = data ); summary(l)$coeff

# Fit model with B
l = lm( C ~ A + B, data = data ); summary(l)$coeff
