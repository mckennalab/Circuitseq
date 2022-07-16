#!/bin/bash

eval "$(conda shell.bash hook)"

export assemblies='<absolute path>/results/assembly_result_with_sample_name/'
export results='<absolute path>/annotated_output/'

conda activate plannotate

for i in $assemblies/* ; do plannotate batch -i $i --html --output $results -s _pLann ; done
