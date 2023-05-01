# Edición de la figura del Anosim
# Solo se requiere editar los nombres de los ejes

"plot.anosim" <-
  function (x, title=NULL, ...) 
  {
    boxplot(x$dis.rank ~ x$class.vec, notch=TRUE, varwidth=TRUE,
            xlab="Grupos", ylab="Rangos de Distancia", cex.lab=1.5,
            ...)
    title(title)
    if (x$permutations) {
      pval <- format.pval(x$signif)
    } else {
      pval <- "not assessed"
    }
    mtext(paste("R = ", round(x$statistic, 3), ", ",
                "P = ", pval ), 3)
    invisible()
  }
