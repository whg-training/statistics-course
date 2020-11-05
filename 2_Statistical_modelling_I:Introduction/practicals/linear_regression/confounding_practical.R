echo <- function( message, ... ) {
    cat( sprintf( message, ... ))
}

variable.names = c( "A", "B", "C" )
A = matrix( 0, nrow = 3, ncol = 3, dimnames = list( variable.names, variable.names ))

# Make a list of all diagrams
# We assume:
# A has effect on C
# B either has effect on A or C or both (otherwise not interesting)
# and 
z01 = c(0,1)
diagrams = expand.grid(
    AA =   0, AB = z01, AC = z01,
    BA = z01, BB =   0, BC = z01,
    CA = z01, CB = z01, CC =   0
)

# Assume causation A->C
diagrams = diagrams[ rowSums( diagrams ) == 2, ]
diagrams = diagrams[ diagrams$AB + diagrams$BA < 2, ]
diagrams = diagrams[ diagrams$AC + diagrams$CA < 2, ]
diagrams = diagrams[ diagrams$BC + diagrams$CB < 2, ]
diagrams = diagrams[ diagrams$BC + diagrams$CB < 2, ]


# Can print a diagram using graphviz:
diagram.to.graphviz <- function( diagram, graph.name = "G", variable.names = c( "A", "B", "C" ) ) {
    str = sprintf( "digraph %s {", graph.name )
    for( v in variable.names ) {
        str = paste( str, sprintf( "%s;", v ), sep = "\n" )
    }
    for( v1 in variable.names ) {
        for( v2 in variable.names ) {
            if( diagram[,sprintf("%s%s",v1,v2)] == 1 ) {
                str = paste( str, sprintf( "%s -> %s", v1, v2 ), sep = "\n" )
            }
        }
    }
    str = paste( str, "}\n", sep = "" )
    return( str )
}

# Requires package DiagrammeR
for( i in 1:nrow( diagrams )) {
    print( grViz( diagram.to.graphviz( diagrams[i,] )) )
    Sys.sleep(0.5) # be nice to browser window
}

diagram.to.matrix <- function( diagram, variable.names = c( "A", "B", "C" ) ) {
    return(
        matrix(
            as.integer(diagram),
            nrow = D, ncol = D,
            byrow = T,
            dimnames = list( variable.names, variable.names )
        )
    )
}

simulate.from.diagram <- function( diagram, N = 1000, effect.size = 0.5, variable.names = c( "A", "B", "C" ) ) {
    # convert the diagram to a matrix
    D = length( variable.names )
    diagram = diagram.to.matrix( diagram )
    echo( "Diagram:\n" )
    print( diagram )

    # We think of each node in the diagram as having an extra parent 'epsilon'.
    # Epsilon is a standard gaussian variable and contributes just enough to make the
    # node have variance 1.
    # The calculation is now based on the total contribution of the epsilons to each node
    # and we keep track of this with the following matrix.
    contributions = matrix(0, nrow = D, ncol = D, dimnames = list( sprintf( "e%s", variable.names ), variable.names ))

    # We compute from the root nodes down
    # and successively expand until contributions to all nodes have been computed
    roots = which( colSums( diagram ) == 0 )
    stopifnot( length(roots) > 0 )
    for( node in roots ) {
        contributions[node,node] = 1 ;
    }
    echo( "Roots are:\n" )
    print( roots )

    # computed nodes contains nodes where all contributions
    # have already been computed.
    computed.nodes = roots ;
    while( length( computed.nodes ) < length( variable.names )) {
        # find nodes that only have arrows from those already fully computed.
        next.nodes = setdiff(
            which( colSums( diagram ) == colSums( diagram[computed.nodes,,drop = FALSE ])),
            computed.nodes
        )
        for( node in next.nodes ) {
            # find all incoming arrows and add contributions
            in.arrows = which( diagram[,node] == 1 )
            for( arrow in in.arrows ) {
                contributions[,node] = contributions[,node] + effect.size * contributions[,arrow]
            }
            contributions[node,node] = sqrt( 1 - sum( contributions[,node]^2 ))
            echo( "Node: %d\n", node )
            print( contributions )
            # node is now fully computed
        }
        computed.nodes = c( computed.nodes, next.nodes )
    }

    epsilon = matrix( rnorm( N * D ), nrow = N, ncol = D, dimnames = list( 1:N, variable.names ))
    
    return(
        list(
            diagram = diagram,
            simulated.data = epsilon %*% contributions,
            covariance = t(contributions) %*% contributions,
            contributions = contributions
        )
    )
}

simulations = list()
for( i in 1:nrow( diagrams )) {
    simulations[[i]] = simulate.from.diagram( diagrams[i,], N = 10000 )
    simulations[[i]]$unadjusted = summary( lm( C ~ A, data = as.data.frame( simulations[[i]]$simulated.data )))$coeff
    simulations[[i]]$adjusted = summary( lm( C ~ A + B, data = as.data.frame( simulations[[i]]$simulated.data )) )$coeff
}
for( i in 1:nrow( diagrams )) {
    echo( "Diagram %d...\n", i )
    print( simulations[[i]]$diagram )
    print( simulations[[i]]$covariance )
    print( simulations[[i]]$unadjusted[,1:2] )
    print( simulations[[i]]$adjusted[1:2,1:2] )
}
