
INNER_DATA=GIarea_byTarget.tsv
OUTER_DATA=L1000metric_byTarget.tsv


pre: include.samples include.features inner.centered.tab outer.centered.tab 

inner.centered.tab : ${INNER_DATA}
	./normalizeForCircles.py $< > $@

outer.centered.tab : ${OUTER_DATA}
	./normalizeForCircles.py $< > $@

include.features:
	cut -f 1 ${INNER_DATA} | grep -v role > $@

include.samples:
	head -n 1 ${INNER_DATA} | transpose.pl > $@

img : 
	mkdir img
	./circlePlot.py -s include.samples -f include.features $@ \
	${INNER_DATA}:color_spec/giarea.color.txt \
	${OUTER_DATA}:color_spec/l1000.color.txt \
	 >& draw_ipl_centered.err &

clean :
	rm -rf img
