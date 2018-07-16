AAACTTAGCTCAAligned.sortedByCoord.out.sorted <- c(0.0,0.0,0.0275229357798,0.0366972477064,0.0688073394495,0.045871559633,0.0596330275229,0.110091743119,0.133027522936,0.146788990826,0.119266055046,0.178899082569,0.284403669725,0.197247706422,0.155963302752,0.220183486239,0.284403669725,0.201834862385,0.252293577982,0.298165137615,0.357798165138,0.302752293578,0.247706422018,0.192660550459,0.183486238532,0.192660550459,0.229357798165,0.288990825688,0.279816513761,0.270642201835,0.206422018349,0.197247706422,0.233944954128,0.288990825688,0.275229357798,0.316513761468,0.316513761468,0.348623853211,0.376146788991,0.371559633028,0.412844036697,0.399082568807,0.394495412844,0.389908256881,0.362385321101,0.353211009174,0.385321100917,0.380733944954,0.385321100917,0.307339449541,0.293577981651,0.339449541284,0.357798165138,0.385321100917,0.330275229358,0.412844036697,0.348623853211,0.385321100917,0.399082568807,0.45871559633,0.463302752294,0.495412844037,0.522935779817,0.48623853211,0.47247706422,0.440366972477,0.5,0.51376146789,0.623853211009,0.577981651376,0.605504587156,0.573394495413,0.550458715596,0.55504587156,0.605504587156,0.628440366972,0.596330275229,0.674311926606,0.683486238532,0.798165137615,0.788990825688,0.862385321101,0.816513761468,0.779816513761,0.733944954128,0.669724770642,0.766055045872,0.908256880734,1.0,0.977064220183,0.93119266055,0.880733944954,0.899082568807,0.839449541284,0.821100917431,0.665137614679,0.545871559633,0.380733944954,0.293577981651,0.123853211009)


pdf("AAACTTAGCTCAAligned.sortedByCoord.out.bam.geneBodyCoverage.curves.pdf")
x=1:100
icolor = colorRampPalette(c("#7fc97f","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f"))(1)
plot(x,AAACTTAGCTCAAligned.sortedByCoord.out.sorted,type='l',xlab="Gene body percentile (5'->3')", ylab="Coverage",lwd=0.8,col=icolor[1])
dev.off()
