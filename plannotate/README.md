We've had a few requests to try to implement [pLannotate](https://github.com/barricklab/pLannotate). Currently we are running into some issues running with singularity resulting in some plasmids causing unexpected errors. 
For the meantime we have a small bash script that should automate the annotation of your assemblies. 

To begin with set up the conda environment for plannotate as described on their github:
```
conda create -n plannotate -c conda-forge -c bioconda plannotate
```

Then modify the `plannotate.sh` bash script in this directory by exporting the path to the annotated results folder using `assemblies` and set a save path of your choosing with `results`.
Then you can simply run the bash script with:
```
bash plannotate.sh
```

I apologise for the quick and dirty fix, but at least this should allow you to use it while I try to figure out the singularity problems. 
