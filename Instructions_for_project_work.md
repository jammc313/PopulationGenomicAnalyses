# Investigating population structure and affinity of unidentified human populations 

## Instructions for project work
For this project you should organize yourselves into groups of 3 or 4. Work together to complete the project, and the end result should be a scientific report (including Abstract, Introduction, Materials & Methods, Results, Discussion, References, Appendices) of approximately 15 pages. You will present your results together on the **15th October**. This is also the deadline for the written report.

### Overview
You have been given genetic data from target populations and reference datasets for comparison. They can be found here:

```
/proj/uppmax2021-2-19/nobackup/HUMAN_POPGEN_PROJECT/IN_DATA
```

You will also find some useful scripts in the `SCRIPTS` folder. Your task is to try to identify the mystery populations you have in the `Unknown.fam`. 

## Basic workflow
The first thing to note is that your unknown samples are "different" somehow. You should look at the genotyping rate (i.e. missingness) and compare this between your unknown and reference datasets. An important part of your task is to figure out what kind of samples you have, and why they are valuable. A further hint is that it can be helpful to look at the sample names in the data files! 

Below is a very rough outline of what you need to do

* Check for related individuals in your reference dataset (and filter them out)
* Unify SNP names across your dataset (by position)
* Filter the **reference** individuals & SNPs (**not** your unknown samples)
* Merge datasets
* PCA
* ADMIXTURE
* Projected PCA 


## Population structure analysis
Before we start here are some basic Unix/Linux commands if you are not used to working in a Unix-style terminal:

### Moving about:
```
    cd – change directory
    pwd – display the name of your current directory
    ls – list names of files in a directory
```
### File/Directory manipulations:
```
    cp – copy a file
    mv – move or rename files or directories
    mkdir – make a directory
    rm – remove files or directories
    rmdir – remove a directory
```
Display file content:
```
	cat - concatenate file contents to the screen or to a file
	less - open file for viewing
```	

### Setup
Create your own working directory within the storage project.
```
cd /proj/uppmax2021-2-19/nobackup/HUMAN_POPGEN_PROJECT/ # move into the Uppmax project for this course
mkdir your_directory_name # e.g. Firstname_Lastname
```
Do not edit the original data files in "IN_DATA" or "SCRIPTS". Copy the contents of "SCRIPTS" into your working directory, and create softlinks to the reference datasets and unknown samples. 

### Filtering and merging population genomic data using PLINK
Link to PLINK site:[https://www.cog-genomics.org/plink2](https://www.cog-genomics.org/plink2)

PLINK is a software for fast and efficient filtering, merging, editing of large SNP datasets, and exports to different outputs directly usable in other programs. It stores huge datasets in compact binary format from which data can be read-in very quickly.
 
Running the program:

The software is already pre-installed on Rackham, you simply have to load it. Working on Rackham type:
```
module load bioinfo-tools 
module load plink/1.90b4.9 
```
Try it:
```
plink
```

#### Basic command structure: 
``` 
plink --filetypeflag filename --commandflag commandspecification --outputfilecommand --out outputfilename
```
For example to do filtering of missing markers at 10% frequency cutoff (reading-in a bed format file called file1, doing the filtering, writing to a ped format file called file2):
```
plink --bfile file1 --geno 0.1 --recode --out file2
```

#### Input formats:
File format summary:
ped format: usual format (.ped and .fam)
.ped contain marker and genotype info and .fam files contain sample info
bed format: binary/compact ped format (.fam .bim and .bed)
(.fam - sample info  .bim - marker info  and .bed - genotype info in binary format
tped format: transposed ped format (.tfam and .tped files)
tfam sample info, .tped marker and genotype info in transposed format
lgen format: long format (see manual, not used that often)

The "Unknown" samples are files containing SNPs from unknown population groups in bed file format. You are going to figure out the ancestry of these population groups during this practical. 

Look at the `.bim` (markers) and `.fam` (sample info) files by typing:
```
less unk1.bim
``` 

do the same for the `.fam` file:
```
less unk1.fam 
```
As mentioned before the `.bim` file store the variant markers present in the data and `.fam` lists the samples. You can try to look at the .bed as well, but this file is in binary format and cannot be visualized as text. If you want to see the specific genotype info you must export the .bed format to a .ped format.

Read in a .bed file dataset and convert to .ped format by typing/pasting in:
```
plink --bfile unk1 --recode --out unk1_ped 
```
#### Look at the info that Plink prints to the screen. How many SNPs are there in the data? How many individuals? How much missing data? Compare this level of missing data to that of a reference dataset. Is there a difference, if so; why do you think that is the case? HINT: Look at the names of the unknown samples!

#### Look at the first few lines of the newly generated .map (sample info) and .ped (marker and genotype info) files using the more command as demonstrated above.

Reading in a .bed/.ped file and converting to .tped:
```
plink --bfile unk1 --recode transpose --out unk1_tped 
plink --file unk1_ped --recode transpose --out unk1_tped 
```
Do you see the difference in the two commands above for reading from a .bed (--bfile) and reading from a .ped (--file) file? Which one takes longer to read-in?

Look at the first few lines of the  .tfam and .tped files by using the `less` command

#### Can you see what symbol is used to encode missing data?

Note: try to always work with .bed files; they are much smaller and take less time to read in. You can see this for yourself by inspecting the difference in file sizes:

```
ls -lh * 
```

Plink can convert to other file formats as well, you can have a look in the manual for the different types of conversions possible.

## Data filtering: Missingness, HWE and MAF
Now we will start to clean our data before further analysis. Take note of how many individuals and markers are being filtered out of your dataset at each stage; this is important information for the report.

Look at the missingness information of each individual and SNP by typing:
```
plink  --bfile unk1 --missing --out test1miss
```
Look at the two generated files by using the `less` command (q to quit):
```
less test1miss.imiss
less test1miss.lmiss
```

The `.imiss` contains the individual missingness and the `.lmiss` the marker missingness. Do you understand the columns of the files? The last three columns are the number of missing, the total, and the fraction of missing markers and individuals for the two files respectively.

We will start our filtering process by filtering for missing data.
First, we filter for marker missingness, we have 1000’s of markers but only 80 individuals, so we would prefer to filter markers rather than individuals.

Paste in the command below to filter out markers with more than 10% missing data
```
plink --bfile unk1 --geno 0.1 --make-bed --out unk2 
```
Look at the screen output, how many SNPs were excluded?

Now we will filter for individual missingness.
Paste in the command below to filter out ind
```
plink --bfile unk2 --mind 0.15 --make-bed --out unk3 
```
Look at the screen output, how many individuals were excluded?

Filtering for minimum allele frequency is not always optimal, especially if you intend later to merge your data with other datasets in which those alleles might be present at higher frequencies. So save this step until your reference dataset is fully merged, then apply a MAF filter:

Filter data for a minimum allele frequency of 1% by pasting in:
```
plink --bfile unk3 --maf 0.01 --make-bed --out unk4 
```
How many SNPs are left?

You can also filter for SNPs out of Hardy-Weinberg equilibrium. Most likely, SNPs out of HWE usually indicates problems with the genotyping. However, to avoid filtering out SNPs that are selected for/against in certain groups (especially when working with case/control data) filtering HWE per group is recommended. After, only exclude the common SNPs that fall out of the HWE in the different groups. However, for reasons of time, you can simply filter the entire dataset for SNPs that aren’t in HWE with a significance value of 0.001
```
plink --bfile unk4 --hwe 0.001 --make-bed --out unk5
```

Look at the screen. How many SNPs were excluded?

If you only what to look at the HWE stats you can do as follows. By doing this command you can also obtain the observed and expected heterozygosities. 
```
plink --bfile unk5 --hardy --out hardy_unk5
```
Look at file hardy_unk5.hwe, see if you understand the output?

There are additional filtering steps that you can go further. The PLINK site manual all the cool commands that you can use to treat your data. Usually, we also filter for related individuals and do a sex-check on the X-chromosome to check for sample mix-ups. Let's filter out related individuals from your reference datasets:
In the `SCRIPTS` folder, you will find a script named `sbatch_KING.sh`. This can be used to run [KING](http://people.virginia.edu/~wc9c/KING/manual.html) and filter out related individuals from your dataset. Take a look inside for instructions on how to run the script. Look at the manual and try and figure out how the software works and think about why this is an important step. After you have run the script, look at the produced output files and figure out how to remove the related individuals from your datases. 

Take note: you should not apply the `maf` filtering until you have your final references dataset! 

=============================================================================


#### Data merging and strand flipping

The `rename_SNP.py` script can be used on the `.bim` files to change the SNP name into the names based on position. This in needed when merging datasets from different SNP chips since the same position can have different names on different chips. When merging the different reference datasets, and then finally your target dataset with the with the reference dataset, you want to make sure to only keep overlapping SNPS.

Usually, when you merge your data with another dataset there are strand issues. The SNPs in the other dataset might be typed on the reverse DNA strand and yours on the forward, or vice versa. Therefore you need to flip the strand to the other orientation for all the SNPs where there is a strand mismatch. One should not flip C/G and A/T SNPs because one cannot distinguish reverse and forward orientation (i.e. C/G becomes G/C unlike other SNPs i.e. G/T which become C/A). Therefore before merging and flipping all A/T and C/G SNPs must be excluded. However, this can be a problem since some of your SNPs in your dataset may be monomorphic when you don't apply the MAF filter. I.E in the bim file they will appear as C 0 (with 0 meaning missing). So you don't know what kind of SNP it is, it can be C G or C T for instance if it is C G it needs to be filtered out but not if it is C T.

Therefore, before merging our data to other datasets it is important to first merge your data with a fake / reference_individual, that you prepare, which is heterozygous at every SNP position. This “fake” reference individual you can easily prepare from the SNP info file you get from the genotyping company or your own genomics processing software (such as Genome Studio from Illumina). You can also prepare it from data downloaded for each SNP from a web-database such as dbSNP. This fake / reference individual has been provided for you in the Schlebusch reference data directory. 

Have you noticed that PLINK sometimes generates a .nof file and in the log file output the following is mentioned:
902 SNPs with no founder genotypes observed
Warning, MAF set to 0 for these SNPs (see --nonfounders)
These are the monomorphic SNPs in your data.

So our first step will be merging with a reference that we prepared from SNP info data beforehand:

Extract your SNPs of interest from the RefInd (remember you filtered out several SNPs already).
```
plink --bfile RefInd1 --extract unk5.bim --make-bed --out RefInd1_ext 
```

Make a list of CG and AT SNPs in your data:
```
sed 's/\t/ /g' RefInd1_ext.bim | grep " C G" >ATCGlist
sed 's/\t/ /g' RefInd1_ext.bim | grep " G C" >>ATCGlist
sed 's/\t/ /g' RefInd1_ext.bim | grep " A T" >>ATCGlist
sed 's/\t/ /g' RefInd1_ext.bim | grep " T A" >>ATCGlist
```

Exclude the CG and AT SNPs from both your reference individuals and target data:
```
plink  --bfile RefInd1_ext --exclude ATCGlist --make-bed --out RefInd1_ext2 
plink  --bfile unk5 --exclude ATCGlist --make-bed --out unk6
```

Merge with RefInd:
```
plink --bfile RefInd1_ext2 --bmerge unk6.bed unk6.bim unk6.fam --make-bed --out MergeRef1  
```

An error is generated because of the strand mismatches. The generated file MergeRef1.missnp
contains the info on the SNPs where there are mismatches - flip the strand of these SNPs in your data:
```
plink --bfile unk6 --flip MergeRef1-merge.missnp --make-bed --out  unk7  
```

Try merging again:
```
plink --bfile RefInd1_ext2 --bmerge unk7.bed unk7.bim unk7.fam --make-bed --out MergeRef2  
```

Now it works. No .nof file is generated which means none of your SNPs are monomorphic anymore (as you would expect!).

Now we will merge our data with a set of reference populations that we get from an already published study. Many of the sites archiving the data provided them in PLINK format as well. For this practical, we selected a few Reference datasets containing different African populations to compare against your Unknown populations.

Look at the .fam files of your reference datasets, do you recognize some of these pops? 

First, extract the SNPs we have in our data from the downloaded RefPops
```
plink --bfile refpops1 --extract MergeRef2.bim --make-bed --out refpops2  
```

Now we will merge our data with the downloaded data:
```
plink --bfile MergeRef2 --bmerge refpops2.bed refpops2.bim refpops2.fam --make-bed --out MergeRefPop1  
```

When you run into strand issues, flip the strands of the refPop datasets to be merged as shown above, and try merging until it works.

Look at your screen output. You will see that the Refpops only contains SNPs that overlap with a small percentage of the SNPs in the Unknown Pops data (~15 000 vs ~95 000). We will now again filter for SNP missingness to exclude all of the extra SNPs in the Unknown data. We want to retain only the overlapping SNPs - if we did not do this, we would find later in our analysis that all the unknown individuals would cluster together separately from the reference dataset. Thatresult would occur due to the non-overlapping sets of SNPs present between the reference and Unknown data, rather than any true genetic dissimilarities!

```
plink --bfile MergeRefPop1 --geno 0.1 --make-bed --out MergeRefPop1fil 
```

How many SNPs are left for your analyses? Keep note of this.

The last thing to do is to extract your fake/Ref_ind from your data:
```
plink --bfile MergeRefPop1fil --remove RefInd1.fam --make-bed --out MergeRefPop2  
```

This is the final files for the next exercise. Rename them:
```
mv MergeRefPop2.bed PopStrucIn1.bed; mv MergeRefPop2.bim PopStrucIn1.bim; mv MergeRefPop2.fam PopStrucIn1.fam 
```

Now you have generated your input files for the next exercise which will deal with population structure analysis. You will look at the population structure of your unknown samples in comparison to the known reference populations.

=============================================================================

#### Population structure inference with PCA and ADMIXTURE

Using ADMIXTURE/PONG and principal component analysis with EIGENSOFT.

Admixture is a similar tool to STRUCTURE but runs much quicker, especially on large datasets.
Admixture runs directly from .bed or .ped files and needs no extra parameter file preparation. You do not specify burnin and repeats, ADMIXTURE exits when it converged on a solution (Delta< minimum value)

First, you have to load the module:

```
module load bioinfo-tools
module load ADMIXTURE/1.3.0

```
A basic ADMIXTURE run looks like this:

```
admixture -s time PopStrucIn1.bed 2
```

This command executes the program with a seed set from system clock time, it gives the input file (remember the extension) and the K value at which to run ADMIXTURE (2 in the previous command).

For ADMIXTURE you also need to run many iterations at each K value, thus a compute cluster and some scripting is useful.

Make a script from the code below to run Admixture for K = 2-6 with 3 iterations at each K value:


```
for i in {2..6};
    do                                                                                      
    for j in {1..3};                                                                                      
        do
        admixture -s time PopStrucIn1.bed ${i} 
        mv PopStrucIn1.${i}.Q PopStrucIn1.${i}.Q.${j};
        mv PopStrucIn1.${i}.P PopStrucIn1.${i}.P.${j};
    done
done
```

Look for a while at the screen output. You will see a short burin, followed by the repeats, and the run stops if delta goes below a minimum value. For K=2 this happens quickly, but the higher Ks take longer. If it takes too long for your liking (it probably will take around 5-10 min) press ctrl-C and copy the output from the folder Data.

Look at the generated output files. What do they contain?

You can quickly look at your admixture output in R by opening R and pasting in the code below. 

```
WD<-getwd()
setwd(WD)

k2_1 <- read.table("PopStrucIn1.2.Q.1")
k2_2 <- read.table("PopStrucIn1.2.Q.2")
k2_3 <- read.table("PopStrucIn1.2.Q.3")
k3_1 <- read.table("PopStrucIn1.3.Q.1")
k3_2 <- read.table("PopStrucIn1.3.Q.2")
k3_3 <- read.table("PopStrucIn1.3.Q.3")
k4_1 <- read.table("PopStrucIn1.4.Q.1")
k4_2 <- read.table("PopStrucIn1.4.Q.2")
k4_3 <- read.table("PopStrucIn1.4.Q.3")
k5_1 <- read.table("PopStrucIn1.5.Q.1")
k5_2 <- read.table("PopStrucIn1.5.Q.2")
k5_3 <- read.table("PopStrucIn1.5.Q.3")
k6_1 <- read.table("PopStrucIn1.6.Q.1")
k6_2 <- read.table("PopStrucIn1.6.Q.2")
k6_3 <- read.table("PopStrucIn1.6.Q.3")

pdf (file ="Admixture_Plot1.pdf", width =25, height = 40, pointsize =10)
par(mfrow=c(15,1))
barplot(t(as.matrix(k2_1)), col=rainbow(2),border=NA)
barplot(t(as.matrix(k2_2)), col=rainbow(2),border=NA)
barplot(t(as.matrix(k2_3)), col=rainbow(2),border=NA)
barplot(t(as.matrix(k3_1)), col=rainbow(3),border=NA)
barplot(t(as.matrix(k3_2)), col=rainbow(3),border=NA)
barplot(t(as.matrix(k3_3)), col=rainbow(3),border=NA)
barplot(t(as.matrix(k4_1)), col=rainbow(4),border=NA)
barplot(t(as.matrix(k4_2)), col=rainbow(4),border=NA)
barplot(t(as.matrix(k4_3)), col=rainbow(4),border=NA)
barplot(t(as.matrix(k5_1)), col=rainbow(5),border=NA)
barplot(t(as.matrix(k5_2)), col=rainbow(5),border=NA)
barplot(t(as.matrix(k5_3)), col=rainbow(5),border=NA)
barplot(t(as.matrix(k6_1)), col=rainbow(6),border=NA)
barplot(t(as.matrix(k6_2)), col=rainbow(6),border=NA)
barplot(t(as.matrix(k6_3)), col=rainbow(6),border=NA)
dev.off()
q()
N

```

This creates the pdf `Admixture_Plot1.pdf`. The bar plots have the individual K cluster assignment for the 3 iterations at K=2-6. The order of individuals is in file “PopStrucIn1.fam”

## PONG 

The method above is a way to quickly check your data, but you have to look at each iteration separately. This makes it hard to get a good overview of the results. We instead combine the different iterations using a software called PONG. 

A typical PONG command looks like this:


```
pong -m your_filemap.txt -i your_ind2pop.txt -n your_pop_order.txt -g
```

To be able to run PONG we thus need to generate three different files.

The first being the filemap. This is the only input that is strictly required to run PONG. IT consists of three columns.
From the PONG manual: 


```
Column 1. The runID, a unique label for the Q matrix (e.g. the string “run5_K7”).

Column 2. The K value for the Q matrix. Each value of K between Kmin and Kmax must
be represented by at least one Q matrix in the filemap; if not, pong will abort.
Column 3. The path to the Q matrix, relative to the location of the filemap. 
```

In order to create what we need we can run the following loop:


```
for i in {2..6};
do
    for j in {1..3};
    do
    echo -e "k${i}_r${j}\t$i\tPopStrucIn1.${i}.Q.${j}" >> unknown_filemap.txt
    done
done
```

The next file we need to create is the ind2pop file. It is just a list of which population each individual belongs to.
We have this information in the `.fam` file so we can just cut out the field we need:

```
cut -f 1 -d " " PopStrucIn1.fam > unknown_ind2pop.txt
``` 


The poporder file is a key between what your populations are called and what "Proper" name you want to show up in your final plot.
For us, it will look like this. *Note that the file needs to be tab-delimited, i.e separated by tabs* 

```
CEU	European
Han	Han_Chinese 
MbutiPygmies	MbutiPygmies
San	San
Unknown1	Unknown1
Unknown11	Unknown11
Unknown3	Unknown3
Unknown5	Unknown5
YRI	Yoruba 

```
So just copy this into a file called `unknown_poporder.txt`

Now we have all the files we need. Time to run PONG.
PONG is available through the module system on Uppmax

```
module load pong 
```

Since we are several people who are going to run PONG at the same time we need to use a different port, otherwise, we will collide with each other. The default port for PONG is 4000. Any other free port will work, like 4001, 2, etc. Make sure you are using a uniq port before proceeding


```
pong -m unknown_filemap.txt -i unknown_ind2pop.txt -n unknown_poporder.txt -g --port YOUR_PORT_NUMBER_HERE
```


When PONG is done it will start hosting a webserver that displays the results at port 4000 by default:  http://localhost:4000. Pong needs to be running for you to look at the results, i.e. if you close it it will not work..

The web server is hosted on the same login node as you were running pong. In case you are unsure of which one that is you can use the command `hostename` to figure it out:

```
hostname
```

To view files interactively you need to have an X11 connection. So when you connect to rackham do:


```
ssh -AY YOUR_USERNAME_HERE@rackham.uppmax.uu.se

```

Make sure that you connect to the same rackham (i.e 1,2,3 etc) as you got from hostename.  
In a new tab (if you didn't put PONG in the background) type:

```
firefox http://localhost:YOUR_PORT_NUMBER
```


#### What do you see? What information does this give you about your unknown populations? 
#### Can you figure out who they are?

Once you are finished with looking at the PONG output you can click and save some output to file and then close the program by hitting `ctrl - c` in the terminal tab running it. 


## Principal component Analysis with Eigensoft

The last population structure method we will look at is Principal Components Analysis with Eigensoft.

You run eigensoft on the .bed format from plink. You only need to modify your .fam file a little for it to work in Eigensoft. The .bed and .map files you use directly as-is. The .fam file you change the extension to .pedind and you substitute the last column (-9 at the moment indicating missing phenotype) with population numbers. When assigning pop numbers do not use 1, 2 or 9. They are reserved for cases, controls and missing data in Eigensoft.

Paste this piece of code in to the terminal:

```
cut -d " " -f1-5 PopStrucIn1.fam >file1a
cut -d " " -f1 PopStrucIn1.fam >file2a
sed "s/Unknown1/51/g" <file2a | sed "s/Unknown3/53/g" | sed "s/Unknown5/55/g" | sed "s/Unknown11/61/g" | sed "s/CEU/81/g" | sed "s/YRI/82/g" | sed "s/Han/83/g" | sed "s/San/84/g" | sed "s/MbutiPygmies/85/g" >file3a
paste file1a file3a >fileComb
sed "s/\t/ /g" fileComb > PopStrucIn1.pedind
rm file1a; rm file2a; rm file3a; rm fileComb
```
It will make a `.pedind` file from your `.fam` file.

Furthermore, you need a parameter file to indicate your parameter options to EIGENSOFT.

Copy the prepared parameter file from the `DATA` directory to your working folder it's called 
`PopStrucIn1.par`

Open the parameter file and look at what is specified in it. At the start is the input/output files. Furthermore, we ask for the info for 10 PCs to be output, qtmode is set to NO to indicate more than one pop, we prune the SNPs based on LD for an r2 value of 0.2. We don't remove any outlying points and we limit the sample size to 20. This is important for PCA where there are groups with very large sample sizes since large sample sizes will distort the PC plot. It is best for PCA that sample sizes are as even as possible.

Run the smartpca package in Eigensoft by typing

```
module load eigensoft

smartpca -p PopStrucIn1.par
```

The outputfiles are .evec and .eval

In the `.evec` file is the main output for the number of PCs that you specified in the .par file. The first row is the Eigenvalues for each of your PCs the rest of the rows list your Pop:Ind specification, the PCs and their loadings and your PopNumber at the end. in the `.eval` is all the eigenvalues that were extracted. To work out the percentage of variation each PC explains, you divide your particular PC eigenvalue by the sum over all the eigenvalues.

We will plot the PCs in R now


Prep for R:

```
sed 1d PopStrucIn1.evec | sed  "s/:/   /g " >   PopStrucIn1.evecm
```


Open `R` and paste the following code to plot your PCs:

```
WD<-getwd()
setwd(WD)
## Define these
evec<- read.table ("PopStrucIn1.evecm")
eval<- read.table ("PopStrucIn1.eval")
namer <- "PopStrucIn1"
nrpc<-10
## Script start
totalev <-sum(eval)
aa <- array(NA,dim=c(nrpc,1))
for (i in 1:nrpc) {
    aa[i,1]<-format(round(((eval[i,1]/totalev)*100),3), nsmall = 3)}
pdf (file =paste(namer, "_PCA1.pdf", sep=""), width =10, height = 15, pointsize =12)
par(mfrow=c(3,2), oma=c(0,0,4,0))
plot (evec[,3], evec[,4],  pch = as.numeric(evec[,1]), cex = 1, col = as.numeric(evec[,1]), xlab=paste("PC1: ", aa[1,1], "%", sep=""), ylab=paste("PC2: ", aa[2,1], "%", sep=""))
legend("topright",  legend = unique(evec[,1]), text.col = "black", cex = 0.75, pch =unique(evec[,1]), col = unique(evec[,1]), xpd = TRUE, bty="n", )
abline(h=0, col="lightgray", lty=5, lwd=0.8); abline(v=0, col="lightgray", lty=5, lwd=0.8)
plot (evec[,5], evec[,6],  pch = as.numeric(evec[,1]), cex = 1, col = as.numeric(evec[,1]), xlab=paste("PC3: ", aa[3,1], "%", sep=""), ylab=paste("PC4: ", aa[4,1], "%", sep=""))
abline(h=0, col="lightgray", lty=5, lwd=0.8); abline(v=0, col="lightgray", lty=5, lwd=0.8)
plot (evec[,7], evec[,8],  pch = as.numeric(evec[,1]), cex = 1, col = as.numeric(evec[,1]), xlab=paste("PC5: ", aa[5,1], "%", sep=""), ylab=paste("PC6: ", aa[6,1], "%", sep=""))
abline(h=0, col="lightgray", lty=5, lwd=0.8); abline(v=0, col="lightgray", lty=5, lwd=0.8)
plot (evec[,9], evec[,10],  pch = as.numeric(evec[,1]), cex = 1, col = as.numeric(evec[,1]), xlab=paste("PC7: ", aa[7,1], "%", sep=""), ylab=paste("PC8: ", aa[8,1], "%", sep=""))
abline(h=0, col="lightgray", lty=5, lwd=0.8); abline(v=0, col="lightgray", lty=5, lwd=0.8)
plot (evec[,11], evec[,12],  pch = as.numeric(evec[,1]), cex = 1, col = as.numeric(evec[,1]), xlab=paste("PC9: ", aa[9,1], "%", sep=""), ylab=paste("PC10: ", aa[10,1], "%", sep=""))
abline(h=0, col="lightgray", lty=5, lwd=0.8); abline(v=0, col="lightgray", lty=5, lwd=0.8)
barplot (as.numeric(aa[,1]), xlab = "PC", ylab = "%Variation explained", axes=TRUE)
title(paste(namer, "PC plot"), outer=TRUE, cex.main = 1)
dev.off()
q()
```

See if you understand the code above. It plots PC1vsPC2 etc and puts labels on the plot.

Look at the output PDF. Do the results of your population PCA correspond to the population structure results you got from the ADMIXTURE plots? How many of the PCs do you think contain useful information. What part of the variation is represented by each of the PCs. Can you see the percentage variation that each PC explains?

We also are going to do a projected PC. We will set-up a run in Eigensoft so that the PCs are calculated based on the Ref pops only and the Unknown individuals are projected on the PCs based on the Ref pops.

To do that we need a slightly modified .par file and an additional file that lists the pops to use as ref pops. Copy these two files from scripts and open them to see how they differ from the .par file used above:
`PopStrucIn1_Proj.par`
RefGroups.groups

Run Eigensoft by typing:

```
smartpca -p PopStrucIn1_Proj.par
```

Prepare the output files for R and run the R script as explained above. Remember the output filenames changed thus adapt the scripts used above accordingly (both the bash script line and the first part of the R script enclosed by the “”)

Look at the output PDF, what has changed? This should give you a further indication as to who your Unknown populations are.

Lastly, we will look at which SNPs contribute to which axes in the PC (SNP weightings of each principal component). This is one way to identify the most informative SNPs that define the structure between certain populations. This is useful if you want to look at population structure but you only want to pick a few best SNPs to type. To generate such a list you just add the snpweightoutname to your par file. We will also not prune for LD so that we do not exclude possible informative SNPs
 
Copy the modified .par file look how it looks and run Eigensoft
PopStrucIn1_snpweight.par

```
./smartpca -p PopStrucIn1_snpweight.par
```

Look at the `PopStrucIn1.snpweight_out` file. The file contains a list of snps and the weights they have on each PC. You can easily select and sort the list to obtain the SNPs with max info for a given PC.

Look at your previous generated PC plot `PopStrucIn1_Proj_PCA1.pdf`
Axis one (PC1) explains ~8% of the variation in your whole Ref_pop dataset. It is the axis that defines the difference between African and non-African populations. To generate a sorted list of the SNPs that would be the best to type to look at the difference between African and non-Africans paste the following command:

```
sed 's/ \+ /\t/g' PopStrucIn1.snpweight_out | sed 's/^\t//g' | cut -f1,3 | sort -r -n -k 3,3  >topSNPsPc1
```

Look at the topSNPsPc1 file that is generated


If you are interested in the frequency of the top SNPs in your data. Copy the two scripts below:

```
Extract_snp_from_bed_make_Barplot
Extract_snp_from_bed_make_Barplot.R
```

Open the main script and paste the rs name of the SNP you are interested in in the place where rs00000000 is at the moment. Save the script and run by typing:

```
./Extract_snp_from_bed_make_Barplot
```

You can run the script multiple times for different top of the list (positive) and bottom of the list (negative) and middle (0) SNPs (each time substitute the particular rs number in the script) and compare the output PDFs.
(The script is not adapted to visualize SNP frequencies for SNPs that contain missing data. So if you get the message
mv: cannot stat `rsxxxx_counts_frequencies.txt': No such file or directory
it means there are missing data for this SNP. Just try another SNP until you get representatives of the top bottom and middle SNPs. Then compare the plots to each other)

What can you see about the frequencies of the top (positive) and bottom of the list (negative) SNPs. What can you say about the SNPs around 0?

If you have time you can look at the second PC as well and visualize the frequencies as described above. To extract and sort info for second PC:

```
sed 's/ \+ /\t/g' PopStrucIn1.snpweight_out | sed 's/^\t//g' | cut -f1,4 | sort -r -n -k 3,3  >topSNPsPc2
 
```
This is the end of section 2. Did you figure out who the unknown populations where?!

##################################################################################

Note that before doing a PCA you should prune your data for SNPs in LD. Think about why this is very important! Go to the admixture manual and search for `prune` and you will find out how to carry out the pruning: https://dalexander.github.io/admixture/admixture-manual.pdf

After that you can run PCA and Admixture, remember that Admixture and Projected PCA both take quite a lot of time! (You will require up to 3 days for each to finish in time).

It is a good idea to run PCA on just your reference dataset as well as the final dataset with both the references and the unkown samples. Look at the difference in projections between the two and think carefully about why any differences may have occurred. 


#### Projected PCA
In `SCRIPTS` there is a script called `prep_run_proj_pca.sh` that can be used to run a projected PCA. (Do you know what a Projected PCA is for? And why it is useful in this case?). The script will require some tweaking and you will need a way to plot it. Contact me when are ready to run it! 


### Further reading
[Schlebusch 2012](https://pubmed.ncbi.nlm.nih.gov/22997136/)

[Gurdasani_2015](https://www.nature.com/articles/nature13997)
 
[Patin 2017](https://science.sciencemag.org/content/356/6337/543)

You will find it worthwhile to read these articles, or at least look at their key findings and figures. This will provide you with some background to past major demographic movements occuring in Africa.

Good luck!

James

 

