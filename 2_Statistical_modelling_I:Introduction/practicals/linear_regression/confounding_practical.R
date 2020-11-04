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
diagrams = expand.grid( AA = 0, AB = c(0, 1), AC = c(0, 1), BA = c(0, 1), BB = 0, BC = c(0, 1), CA = 0, CB = 0, CC = 0 )

# Assume causation A->C
diagrams = diagrams[ diagrams$AC == 1, ]
# Allow no reverse causation from C to A
diagrams = diagrams[ diagrams$CA == 0 & ( diagrams$CB + diagrams$BA ) < 2, ]
# Also don't allow causal loops
diagrams = diagrams[ diagrams$BC + diagrams$CB < 2, ]
diagrams = diagrams[ diagrams$AB + diagrams$BA < 2, ]

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

for( i in 1:nrow( diagrams )) {
    print( grViz( diagram.to.graphviz( diagrams[i,] )) )
    Sys.sleep(0.5) # be nice to browser window
}



simulate.from.diagram <- function( diagram, N = 1000, effect.size = 0.5, variable.names = c( "A", "B", "C" ) ) {
    # convert the diagram to a matrix
    D = length( variable.names )
    diagram = matrix( as.integer(diagram), nrow = D, ncol = D, byrow = T, dimnames = list( variable.names, variable.names ))

    # we'll simulate N samples
    # We return simulated data, plus exact covariance matrix
    # We have to use an environment here because R has copy-on-write semantics
    env = new.env()
    env$simulated.data = matrix( 0, nrow = N, ncol = length( variable.names ), dimnames = list( 1:N, variable.names ))
    env$covariance = matrix( 0, nrow = D, ncol = D, dimnames = list( variable.names, variable.names ))

    # compute node fills up the node to variance 1
    # traces this through to descendants
    compute.node <- function( node, effect.size, result ) {
        # recursively update simulated data
        update.descendants <- function( node, contribution, effect.size ) {
            children = which( diagram[node,] == 1 )
            for( child in children ) {
                child.contribution = effect.size * contribution
                env$simulated.data[,child] = env$simulated.data[,child] + child.contribution
                update.descendants( child, child.contribution, effect.size )
            }
        }
        
        # recursively compute covariance matrix
        compute.covariance.from.source <- function( source, node, node.variance, effect.size ) {
            children = which( diagram[node,] == 1 )
            for( child in children ) {
                child.variance = sqrt(effect.size) * node.variance
                echo( "Updating node %d from source %d, child variance = %.2f\n", child, source, child.variance )
                env$covariance[source,child] = env$covariance[source,child] + child.variance
                compute.covariance.from.source( source, child, child.variance, effect.size )
            }
        }

        remaining.variance = 1 - env$covariance[node,node]
        epsilon = rnorm( N, sd = sqrt( remaining.variance ))
        env$simulated.data[,node] = env$simulated.data[,node] + epsilon
        update.descendants( node, epsilon, effect.size )
        env$covariance[node,node] = 1
        compute.covariance.from.source( node, node, 1, effect.size )
    }

    compute.nodes <- function( nodes, effect.size ) {
        for( node in nodes ) {
            compute.node( node, effect.size )
        }
        return( nodes )
    }


    # Every node in our graph is either a root (no arrows in)
    # a first-level child (direct arrow from root)
    # or a 2nd level child (direct arrow from 1st level child)
    # find nodes not caused by anything else
    roots = which( colSums( diagram ) == 0 )
    stopifnot( length(roots) > 0 )
    echo( "Roots are:\n" )
    print( roots )

    computed.nodes = compute.nodes( roots, effect.size )
    while( length( computed.nodes ) < length( variable.names )) {
        # find nodes that only have arrows from those already fully computed.
        next.nodes = setdiff( which( colSums( diagram ) == colSums( diagram[computed.nodes, ])), computed.nodes )
        echo( "Computed nodes:\n" )
        print( computed.nodes )
        echo( "Next nodes:\n" )
        print( next.nodes )
        computed.nodes = c( computed.nodes, compute.nodes( next.nodes, effect.size ))
        print( env$covariance )
    }
    result = list(
        covariance = env$covariance,
        simulated.data = env$simulated.data
    )
    rm( env )
    return( result )
}
