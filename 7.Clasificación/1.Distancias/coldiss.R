# coldiss()
# Color plots of a Matriz de Disimilitud, without and with ordering
#
# License: GPL-2 
# Author: Francois Gillet, August 2009
#

"coldiss" <- function(D, nc = 4, byrank = TRUE, diag = FALSE)
{
	require(gclus)

	if (max(D)>1) D <- D/max(D)

	if (byrank) {
		spe.color = dmat.color(1-D, cm.colors(nc))
	}
	else {
		spe.color = dmat.color(1-D, byrank=FALSE, cm.colors(nc))
	}

	spe.o = order.single(1-D)
	speo.color = spe.color[spe.o,spe.o]
	
	op = par(mfrow=c(1,2), pty="s")

	if (diag) {
		plotcolors(spe.color, rlabels=attributes(D)$Labels, 
			main="Matriz de Disimilitud", 
			dlabels=attributes(D)$Labels)
		plotcolors(speo.color, rlabels=attributes(D)$Labels[spe.o], 
			main="Matriz de Disimilitud Ordenada", 
			dlabels=attributes(D)$Labels[spe.o])
	}
	else {
		plotcolors(spe.color, rlabels=attributes(D)$Labels, 
			main="Matriz de Disimilitud")
		plotcolors(speo.color, rlabels=attributes(D)$Labels[spe.o], 
			main="Matriz de Disimilitud Ordenada")
	}

	par(op)
}

# Usage:
# coldiss(D = dissimilarity.matrix, nc = 4, byrank = TRUE, diag = FALSE)
# If D is not a Matriz de Disimilitud (max > 1), then D is divided by max(D)

# Example:
# coldiss(spe.dj, nc=9, byrank=F, diag=T)

# byrank= TRUE		equal-sized categories
# byrank= FALSE		equal-length intervals
