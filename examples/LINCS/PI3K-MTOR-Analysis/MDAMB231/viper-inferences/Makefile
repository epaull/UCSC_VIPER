

regulators_compare:
	../../scripts/splitSubtypes.py  ../../shared_data/LJP56_Drugs_annotated.tsv viperInferences.txt
	../../scripts/compareRegulators.R --class_1 samples_AktPI3K.txt --class_2 samples_MAPK.txt --regulator FOXO1
	../../scripts/compareRegulators.R --class_1 samples_AktPI3K.txt --class_2 samples_MAPK.txt --regulator FOXO3
	../../scripts/compareRegulators.R --class_1 samples_AktPI3K.txt --class_2 samples_MAPK.txt --regulator FOXO4
