# Barcoding 各软件的优缺点

[TOC]

**First of all!!!**
```python
'重要，对于barcoding研究使用，尽量做到同种序列有多个，否则影响很大'
大部分软件对于未鉴定正确的结果，无论如何都会聚类到一种上，并不会给出'未鉴定'的结果，所有结果可以作为鉴定参考
```


## SpeciesIdentifiter
- 结果查看**列并没有对齐，需要手动调整下，具体位置见下表**

| Query                      | 以第一列为参照，同种，种内 | intraspeciesDist | seqLength | 以第一列为参照，不同种，种间 | interspeciesDist | seqLength | Match                      | Identification                                |
| -------------------------- | -------------------------- | ---------------- | --------- | ---------------------------- | ---------------- | --------- | -------------------------- | --------------------------------------------- |
| TBC03,Tachina nupta        | TBC86,Tachina nupta        | 0.85             | 590       | TBC56,Peleteria iavana       | 8.42             | 616       | TBC86,Tachina nupta        | Successful match at 0.85% (within threshold)  |
| TBC09,Exorista hyalipennis | TBC18,Exorista hyalipennis | 3.09             | 597       | TBC106,Bessa parallela       | 9.06             | 621       | TBC18,Exorista hyalipennis | Successful match at 3.09% (outside threshold) |
| TBC106,Bessa parallela     | TBC84,Bessa parallela      | 0                | 621       | TBC09,Exorista hyalipennis   | 9.06             | 621       | TBC84,Bessa parallela      | Successful match at 0.0% (within threshold)   |
| TBC110,Blepharipa zebina   | TBC129,Blepharipa zebina   | 4.52             | 621       | TBC71,Blepharipa latigena    | 6.23             | 621       | TBC129,Blepharipa zebina   | Successful match at 4.52% (outside threshold) |
| TBC115,Trixa longipennis   | TBC53,Trixa longipennis    | 1.96             | 621       | TBC86,Tachina nupta          | 14.14            | 590       | TBC53,Trixa longipennis    | Successful match at 1.96% (within threshold)  |




## BarcodingR

> 在读取文件的时候，需要将序列文件读取为DNAbin格式，且只能通过ape包的`read.dna("filename", format = "fasta")`读取，或通过**adegenet**包的`fasta2DNAbin("filename")`来读取，其它方式读取的序列文件不能被**BarcodingR**所识别。
>
> 只能是ref数据中有的种才能被成功鉴定，对于没有参考的数据无法鉴定
>
> **格式强调** fasta格式名称中不能含有","，BarcodingR中除外，其名称为：`>TBC000,Tachinidae_Tachina_nupta`

```R
#library(ape)
#library(adegenet)
#library(phyloch)
#library(BarcodingR)
# for large fasta file, read.fas() package:phyloch is suggest

#ref <- read.dna("ref.fasta", format = "fasta") # read fasta file from your work directory
#que <- fasta2DNAbin("que.fasta") #not recommand, slow

list <- sample.ref(read.dna(file.choose(),format = "fasta"),sample.level = "species")
ref <- list$ref.selected
que <- list$ref.left
bp <- barcoding.spe.identify(ref, que, method = "bpNewTraining") #method = "fuzzyId" or "Bayesian"
bp1 <- barcoding.spe.identify(ref, que, method = "fuzzyId")
bp2 <- barcoding.spe.identify(ref, que, method = "Bayesian")
# compare three methods
# bp <- identified results by bpNewTraining
# bp1 <- identified results by fuzzyId
# bp2 <- identified results by Bayesian
########################################
que.IDs <- gsub(",.+","",rownames(que)) # extract seqs names
bpid <- bp$output_identified$spe.Identified # get values from bpNewTraining IDENTIFICATION
fuzzyid <- bp1$output_identified$spe.Identified # get values from fuzzyId
bayesid <- bp2$output_identified$output_identified # get values from Bayesian
identifications <- data.frame(queIDs = que.IDs, bpid = bpid, fuzzyid = fuzzyid, Bayesianid = bayesid) #create dataframe from three identification \ bpid maybe pid
ccs <- consensus.identify(identifications) # consensus test
save.ids(outfile = "bpNewTraing.txt", bp)
save.ids(outfile = "fuzzyId.txt", bp1)
save.ids(outfile = "Bayesian.txt", bp2)
write.csv(ccs, file = "ccs.csv")
```


**The BarcodingR input fasta format like:**


```
>TBC06,Nilea_hortulana 
GGAGCTTGATCAGGAATAATTGGTACTTCTTTAAGTATTTTAATCCGAACTGAATTAGGACATCCTGGATCTTTAATTGGAGATGACCAAATCTATAATGTAATTGTTACAGCTCATGCTTTCATCATAATTTTTTTTATAGTTATACCAATTATAATTGGAGGGTTTGGAAATTGATTAATTCCTTTAATATTAGGAGCTCCTGATATAGCTTTCCCACGAATAAATAATATAAGTTTTTGATTACTTCCCCCTGCTTTAACACTTTTGTTGGTAAGTAGAATAGTAGAAAATGGAGCTGGTACTGGATGAACAGTTTACCCACCCTTATCATCTATTATTGCACATGGAGGAGCTTCTGTAGATTTAGCTATTTTTTCCCTTCATTTAGCTGGAATTTCATCTATTTTAGGAGCTGTAAATTTTATTACAACTGTAATTAATATACGATCAACCGGAATTACATTTGATCGAATACCTTTATTTGTTTGATCAGTTGTTATTACAGCCCTATTACTTTTATTATCTTTACCTGTTTTAGCCGGAGCTATTACTATATTATTAACAGATCGAAATTTAAATACTTCCTTTTTTGACCCCGCAGGAGGAGGAGACCCAATT------------------
```



## ABGD ([Automatic Barcode Gap Discovery](http://www.ncbi.nlm.nih.gov/pubmed/21883587))

`Note：对于ABGD参数定义的relative length仍需进一步理解，Pmin和Pmax代表种内遗传距离范围`



## GMYC ([general mixed Yule-coalescent model](https://bmcevolbiol.biomedcentral.com/articles/10.1186/1471-2148-12-196))

> Pons, J., et al., Sequence-based species delimitation for the DNA taxonomy of undescribed insects. Syst Biol, 2006. 55.

使用Beast软件构建的超度量树作为输入

**Pathd8**用来校准非超度量树为超度量树，最好以C源码自己编译运行。 [下载地址](http://www2.math.su.se/PATHd8/)  

**R**版本的超度量树校准 *adjust_tip_lengths.r* 。见下方代码块：

`Type **cc PATHd8.c -O3 -lm -o PATHd8 **to compile. Notice that O in -O3 is a letter, not a number.`

GMYC网页版和R相同 [GMYC web sever](http://species.h-its.org/gmyc/)    [GMYC介绍](http://barralab.bio.ic.ac.uk/downloads.html)

```R
#这是在R语言中执行GMYC的代码
#install.packages("splits", repos="http://R-Forge.R-project.org") #install r package
#install.packages("ape")
library(ape)
library(splits)
tree<-read.nexus(file.choose())
test<-gmyc(tree, method="single")
#test<-gmyc(tree, method="multiple")
summary(test)
plot(test)
spec_list <- spec.list(test)
write.table(spec_list, file="GMYCspe_haplotype.txt")
```



### adjust_tip_lengths.r ----get ultrametric tree by R script

```R
## an R script come from Chufei Tang for transfer tree branch length to ultrametric tree
library(paleotree)
library(adephylo)
tree<-read.nexus('parrots_mean_equal_vartim0.1.tre')
lengths = distRoot(tree, tree$tip.label)
desired_length = max(lengths)
terms <- tree$edge[, 2] <= Ntip(tree)
terminal.edges <- tree$edge.length[terms]
names(terminal.edges) <- tree$tip.label[tree$edge[terms, 2]]
for (i in 1:length(tree$tip.label)){
    terminal.edges[tree$tip.label[i]] = terminal.edges[tree$tip.label[i]] + (desired_length-lengths[tree$tip.label[i]])
}
tree$edge.length[terms] = terminal.edges
write.tree(tree,file="parrot_um.tre")
lengths2 = distRoot(tree, tree$tip.label)
```







## PTP&bPTP ([Poisson Tree Processes](https://sco.h-its.org/exelixis/web/software/PTP/index.html))

Code have syntax wrong! *zhangjiajie/PTP: PTP model for species delimitation* [这里是网址](https://github.com/zhangjiajie/PTP)

上面的代码运行有错误，改为网页版 [网址](http://species.h-its.org/ptp/)

> 参考：蔡彦朋博士论文，导师：周红章

> **原理** 它假设每个单独的基因突变导致新物种形成的概率很小，所以相同等位基因在不同物种中的差异应该显著大于种间的差异。该模型整合了两类相互独立的泊松过程（一类用来描述物种形成，另一类描述种内的种群分化），采用最大似然法，通过检索种间分支和种内分支在时间上的转折点来划分物种。bPTP是PTP的升级版，新模型采用贝叶斯和蒙特卡罗算法取代了原来的最大似然法，可以检索更大的解空间。理论上，PTP模型的检索结果是bPTP结果的子集（Zhang *et al*., 2013; Dumas *et al*., 2015).  [引自蔡彦朋，2016博士论文]



## BLOG (2013) [软件下载地址](http://dmb.iasi.cnr.it/blog-downloads.php)  [文章下载](http://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12073/full)

**Description:**

BLOG is a application devoted to the automatic classification of animal species through the analysis of a small portion of mitochondrial DNA, [DNA Barcode](http://www.barcoding.si.edu/). The application is described in detail in Bertolazzi, Felici, Weitschek [Learning to classify species with barcodes](http://www.biomedcentral.com/1471-2105/10/S14/S7). To run BLOG, you must provide a training set of barcodes (650 sites with A, C, G, T) and the species to which each barcode belongs. Then, logic formulas are extracted from the training data and a test set of one of more barcode is classified according to these formulas. 
The input files are standard FASTA barcode sequences, described [here](http://dmb.iasi.cnr.it/faq.php#inputBlog).
The parameters needed for the correct running of BLOG are:

1. Train File - A FASTA file to train BLOG
2. Test File - A FASTA file containing query sequences that require identification fasta file; It will create a train file of almost the 80% of the FASTA sequences, and a train file with the remaining sequences.
3. BLOG fasta format, and  "|" is important

```
>TBC07|Linnaemya_omega
GGAGCTTGATCAGGAATAATTGGAACTTCATTAAGTATCTTAATCCGAGCTGAATTAGGTCATCCAGGTGCATTAATTGGTGATGACCAAATTTATAATGTAATTGTTACAGCCCATGCTTTTGTTATGATTTTTTTTATAGTTATACCAATTATAATTGGAGGGTTCGGAAATTGATTAGTTCCTTTAATACTAGGAGCTCCAGATATAGCCTTTCCACGAATAAATAATATAAGATTTTGACTTTTACCTCCAGCATTAACTTTATTATTAGTAAGTAGCATAGTAGAAAACGGAGCTGGTACAGGATGAACTGTTTATCCACCTTTATCTTCTAACATCGCTCATGGAGGAGCATCTGTTGATTTAGCTATTTTTAGTTTACACTTAGCTGGAATTTCTTCAATTTTAGGAGCCGTAAATTTTATTACTACAGTAATTAATATACGATCTACAGGTATTACTTTTGACCGAATACCTTTATTTGTTTGATCTGTAGTAATCACAGCTTTATTGCTATTATTATCTTTACCAGTATTAGCAGGAGCTATTACAATATTATTAACTGATCGAAATTTAAATACTTCATTCTTTGATCCAGCTGGAGGAGGAGATCCAATT
```



## !*CAOS基于特征的方法* !

2018-3-6国内打不开网址！[下载CAOS](http://sarkarlab.mbl.edu/CAOS) [另一个稍有升级的CAOS](http://www.genomecurator.org/CAOS/P-Gnome/PGnomeindex.html)

> 类似的还有BLOG软件，但经常报错


