load.genes <- function( condense = TRUE ) {
    # Load genes
    gene <- read.delim(
        "resources/refGene_chr19.txt",
        header=TRUE,
        as.is=TRUE
    );
	gene <- gene[ -grep( "hap", gene$chrom ), ]
	
	if( condense ) {
	    gene <- gene[order(gene$txEnd - gene$txStart,decreasing=TRUE),];  #Get just longest transcript
	    gene <- gene[ !duplicated( gene$name2 ), ];
	}
    gene <- gene[ !is.na(gene$txStart), ];
	
    chromosome =  gsub( "^chr", "", gene$chrom )
    w1 = which( nchar( chromosome ) == 1 )
    chromosome[ w1 ] = sprintf( "0%s", chromosome[w1] )
    gene$chromosome = chromosome
    return( gene ) ;
}

plot.genes <- function( genes, region, height_in_inches = 1, exons = get.exons( genes ), vertical = FALSE, ... ) {
	get.exons <- function( genes ) {
		result = data.frame()
		for( i in 1:nrow( genes )) {
			gene = genes[i,]
			if( gene$exonCount > 0 ) {
				result = rbind(
					result,
					data.frame(
						name2 = gene$name2,
						txStart = gene$txStart,
						txEnd = gene$txEnd,
						exonStart = as.integer( strsplit( gene$exonStarts, split = "," )[[1]] ),
						exonEnd = as.integer( strsplit( gene$exonEnds, split = "," )[[1]] )
					)
				)
			}
		}
		return( result )
	}

	w = which( genes$txEnd >= region[1] & genes$txStart <= region[2] )
	print(w)
    if( length(w) > 0 ) {
		genes = genes[w,]
		stopifnot( nrow( genes ) > 0 )
        genes = genes[ order( genes$txStart ),, drop = FALSE ]
        genes$y = NA ;
        genes$y[1] = 1 ;
        if( nrow( genes ) > 1 ) {
            spacer = ( region[2] - region[1] ) / 10 ;
            maxes = ( genes[1,]$txEnd + spacer )
            for( i in 2:nrow( genes )) {
                for( l in 1:length( maxes )) {
                    if( genes$txStart[i] >= maxes[l] ) {
                        genes$y[i] = l ;
                        maxes[l] = genes$txEnd[i] + spacer ;
                        break ;
                    }
                }
                if( is.na( genes$y[i] )) {
                    maxes = c( maxes, genes$txEnd[i] + spacer )
                    genes$y[i] = length( maxes ) ;
                }
            }
        }
        
        exons$y = genes$y[ match( exons$name2, genes$name2 )]
		print( exons )
        height_of_genes = height_in_inches / max( genes$y )
        
        strands = c("+","-" ) ;
        level = c( 80, 70 ) ;
		if( vertical ) {
			plot( region[1], 0, pch = '', ylim = region, xlim = c( 0, max( 2, max(genes$y)+1 ) ), xlab = '', ylab = '', xaxt = 'n', ... ) ;
		} else {
			plot( 0, region[1], pch = '', xlim = region, ylim = c( 0, max( 2, max(genes$y)+0.5 ) ), xlab = 'Genes', ylab = '', yaxt = 'n', ... ) ;
		}

        arrow.sep = ( region[2] - region[1] ) / 200 ;
		arrow.voffset = 0
		
        relative.lengths = ( genes$txEnd - genes$txStart ) / ( region[2] - region[1] )

        genes$mark1 = ( 0.25 * genes$txStart + 0.75 * genes$txEnd ) ;
        genes$mark2 = ( 0.5 * genes$txStart + 0.5 * genes$txEnd ) ;
        genes$mark3 = ( 0.75 * genes$txStart + 0.25 * genes$txEnd ) ;
        genes$sign = 1 ;
        genes$sign[ which( genes$strand == '-' ) ] = -1 ;

		if( vertical ) {
	        segments(
	            y0 = genes$txStart, y1 = genes$txEnd,
	            x0 = genes$y, x1 = genes$y,
	            col = "blue"
	        )
	        segments(
	            y0 = genes$txStart, y1 = genes$txStart,
	            x0 = genes$y + 0.5, x1 = genes$y - 0.5,
	            col = "blue"
	        )
	        segments(
	            y0 = genes$txEnd, y1 = genes$txEnd,
	            x0 = genes$y + 0.5, x1 = genes$y - 0.5,
	            col = "blue"
	        )
		} else {
	        segments(
	            x0 = genes$txStart, x1 = genes$txEnd,
	            y0 = genes$y, y1 = genes$y,
	            col = "blue"
	        )
	        segments(
	            x0 = genes$txStart, x1 = genes$txStart,
	            y0 = genes$y + 0.5, y1 = genes$y - 0.5,
	            col = "blue"
	        )
	        segments(
	            x0 = genes$txEnd, x1 = genes$txEnd,
	            y0 = genes$y + 0.5, y1 = genes$y - 0.5,
	            col = "blue"
	        )
		}
		
        wBigEnough = which( ( genes$txEnd - genes$txStart ) > arrow.sep * 2 ) ;
		if( length(wBigEnough) > 0 ) {
			arrows = data.frame()
			for( i in wBigEnough ) {
				arrows = rbind(
					arrows,
					data.frame(
						name2 = genes$name2[i],
						y = genes$y[i],
						x = seq( from = genes$txStart[i] + arrow.sep, to = genes$txEnd[i], by = arrow.sep ),
						sign = genes$sign[i]
					)
				)
			}
			
			if( vertical ) {
				segments(
					y0 = arrows$x, y1 = arrows$x + arrows$sign * 0.5 * arrow.sep,
					x0 = arrows$y + 0.2 + arrow.voffset, x1 = arrows$y + arrow.voffset,
					col = rgb( 0.5, 0.5, 0.5, 0.7 )
				)

				segments(
					y0 = arrows$x, y1 = arrows$x + arrows$sign * 0.5 * arrow.sep,
					x0 = arrows$y - 0.2 + arrow.voffset, x1 = arrows$y + arrow.voffset,
					col = rgb( 0.5, 0.5, 0.5, 0.7 )
				)

				wExon = which( exons$name2 %in% genes$name2[wBigEnough] )
				if( length(wExon) > 0 ) {
					rect(
						ybottom = exons$exonStart[wExon],
						ytop = exons$exonEnd[wExon],
						xleft = exons$y[wExon] - 0.25,
						xright = exons$y[wExon] + 0.25,
						border = NA,
						col = "blue"
					)
				}
			} else {
				segments(
					x0 = arrows$x, x1 = arrows$x + arrows$sign * 0.5 * arrow.sep,
					y0 = arrows$y + 0.2 + arrow.voffset, y1 = arrows$y + arrow.voffset,
					col = rgb( 0.5, 0.5, 0.5, 0.7 )
				)

				segments(
					x0 = arrows$x, x1 = arrows$x + arrows$sign * 0.5 * arrow.sep,
					y0 = arrows$y - 0.2 + arrow.voffset, y1 = arrows$y + arrow.voffset,
					col = rgb( 0.5, 0.5, 0.5, 0.7 )
				)

				wExon = which( exons$name2 %in% genes$name2[wBigEnough] )
				if( length(wExon) > 0 ) {
					rect(
						xleft = exons$exonStart[wExon],
						xright = exons$exonEnd[wExon],
						ybottom = exons$y[wExon] - 0.25,
						ytop = exons$y[wExon] + 0.25,
						border = NA,
						col = "blue"
					)
				}
			}
		}
		# We fit about 7 gene names per inch at cex = 1.  After that we need to start scaling.
		text.cex = min( 1, ( 7 * height_in_inches ) / max( genes$y ) )
		if( vertical ) {
			text( genes$y, genes$txEnd, label = genes$name2, adj = c( -0.1, 0.5 ), cex = text.cex, srt=90 )
		} else {
			text( genes$txEnd, genes$y, label = genes$name2, adj = c( -0.1, 0.5 ), cex = text.cex )
		}
    } else {
        plot.new()
        plot( region[1], 0, pch = '', xlim = region, ylim = c( 0, (max(genes$y)+1) ), xlab = '', ylab = 'genes' ) ;
    }
}
