# Investigating population structure and affinity of unidentified human populations 

## Instructions for project work
For this project you should organize yourselves into groups of 3 or 4. Work together to complete the project, and the end result should be a scientific report (including Abstract, Introduction, Materials & Methods, Results, Discussion, References, Appendices) of approximately 15 pages. You will present your results together on the **15th October** (which is also the deadline for the written report).

### DATA
You have been given some reference and target populations. They can be found here:

```
/proj/uppmax2021-2-19/nobackup/HUMAN_POPGEN_PROJECT/IN_DATA
```

You will also find some useful scripts in the `SCRIPTS` folder.

## Your task
Ultimately you need to identify what populations you have in the `Unknown.fam`. To achieve this you have been given some reference datasets that should be useful for comparison.



## Workflow
Your unknown samples are "different" in a way, a hint is to look at the genotyping rate - i.e. missingness. Part of your first task is to figure out what types of samples you have, and why they are valuable even though they look the way they do. It can also be helpful to look at their names in the data files. You should not apply the filters to your unknown samples! 

Below is a rough outline of what you need to do

* Check for related individuals in your reference dataset (and filter them out)
* Unify SNP names across your dataset (by position)
* Filter the reference individuals & SNPs
* PCA
* ADMIXTURE
* Projected PCA 


## More thorough instructions

#### Filtering out related individuals
In the `SCRIPTS` folder, there is a script called `sbatch_KING.sh` that can be used to run [KING](http://people.virginia.edu/~wc9c/KING/manual.html). Have a look inside it for instructions on how to run the script. Look at the manual and try and figure out how the software works. After you have run the script, have a look at the produced output files and figure out how to remove the related individuals.

#### Filtering
For filtering procedures, I suggest you take a look at the instructions of this [popgen lab](https://github.com/Hammarn/Populationsgenomik/blob/master/1BG508.md), as well as the other scripts supplied to you. Take note: you should not apply the `maf` filtering until you have your final references dataset! 

```
--mind 0.01
--geno 0.01
--hwe 0.0000001
--maf 0.01
```
Your Unknown samples should not go through any filtering.

#### Merging
The `rename_SNP.py` script can be used on the `.bim` files to change the SNP name into the names based on position. This in needed when merging datasets from different SNP chips since the same position can have different names on different chips. When merging your reference dataset you want to make sure to only keep overlapping SNPS.

When you have merged together your reference dataset and finished filtering it, you can merge it together with your `unkown` data.



#### Running PCA and ADMIXTURE
There are some guidelines for running PCA at the popgen lab link above, and in the SCRIPTS folder. Note that before doing a PCA you should prune your data for SNPs in LD. Go to the admixture manual and search for `prune` and you will find out how to carry out the pruning:
https://dalexander.github.io/admixture/admixture-manual.pdf

After that you can run PCA and Admixture, remember that Admixture takes quite a lot of time!

It's probably a good idea to run PCA on just your reference dataset as well as the final dataset with both the references and the unkown samples.


#### Projected PCA
In `SCRIPTS` there is a script called `prep_run_proj_pca.sh` that can be used to run a projected PCA. 
It will require some tweaking and you will need a way to plot it. Contact me when are ready to run it! 


#### Admixture dating

I highly suggest that you use DATES for admixture dating (but if you prefer to use an alternative, go ahead!):
[https://github.com/priyamoorjani/DATES](https://github.com/priyamoorjani/DATES). It is available through the uppmax module systems as `DATES`. It follows the EIGENSOFT package suites format and syntax.



### Further reading
[Schlebusch 2012](https://pubmed.ncbi.nlm.nih.gov/22997136/)

[Gurdasani_2015](https://www.nature.com/articles/nature13997)
 
[Patin 2017](https://science.sciencemag.org/content/356/6337/543)

You might find it worthwhile to read these articles, or at the very least, look at their key findings and figures.

Good luck!

James

 

