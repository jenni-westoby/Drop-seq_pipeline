AAAACCGGGCATAligned.sortedByCoord.out.sorted <- c(0.0,0.0714285714286,0.166666666667,0.166666666667,0.404761904762,0.714285714286,0.690476190476,0.404761904762,0.547619047619,0.5,0.452380952381,0.619047619048,0.642857142857,0.309523809524,0.404761904762,0.333333333333,0.595238095238,0.714285714286,0.690476190476,0.642857142857,0.452380952381,0.404761904762,0.47619047619,0.595238095238,0.642857142857,0.5,0.547619047619,0.428571428571,0.309523809524,0.452380952381,0.452380952381,0.333333333333,0.166666666667,0.190476190476,0.0714285714286,0.261904761905,0.261904761905,0.357142857143,0.238095238095,0.333333333333,0.333333333333,0.380952380952,0.404761904762,0.571428571429,0.333333333333,0.261904761905,0.404761904762,0.309523809524,0.428571428571,0.166666666667,0.190476190476,0.333333333333,0.190476190476,0.238095238095,0.309523809524,0.357142857143,0.285714285714,0.309523809524,0.238095238095,0.47619047619,0.404761904762,0.547619047619,0.333333333333,0.285714285714,0.238095238095,0.190476190476,0.261904761905,0.309523809524,0.404761904762,0.547619047619,0.428571428571,0.571428571429,0.571428571429,0.452380952381,0.452380952381,0.52380952381,0.428571428571,0.404761904762,0.404761904762,0.309523809524,0.261904761905,0.357142857143,0.452380952381,0.571428571429,0.571428571429,0.571428571429,0.595238095238,0.404761904762,0.285714285714,0.285714285714,0.238095238095,0.404761904762,0.738095238095,0.809523809524,1.0,0.642857142857,0.547619047619,0.642857142857,0.238095238095,0.0952380952381)


pdf("AAAACCGGGCATAligned.sortedByCoord.out.bam.geneBodyCoverage.curves.pdf")
x=1:100
icolor = colorRampPalette(c("#7fc97f","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f"))(1)
plot(x,AAAACCGGGCATAligned.sortedByCoord.out.sorted,type='l',xlab="Gene body percentile (5'->3')", ylab="Coverage",lwd=0.8,col=icolor[1])
dev.off()
