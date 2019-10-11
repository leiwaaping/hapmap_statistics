# hapmap_statistics
this script is to convet HAPMAP format to pca-valued foramt and calculate specific frequency.

heatmap_calculationt
intput：3-testcase.txt 3-testcontrol.txt   
output: two files ,test_matrix.txt is what we want.  
command:  
```python2 hapmap_calculation.py 3-testcase.txt 3-testcontrol.txt test.txt test_matrix.txt```

processing:
```
$ head 3-testcase.txt
rs# alleles NA18532 NA18605 NA18542 NA18550 NA18608 NA18564 NA18571 NA18620 NA18623 NA18633 NA18637 NA18526 NA18524 NA18558 NA18562 NA18537 NA18603 NA18545 NA18572 NA18547 NA18609 NA18552 NA18555 NA18566 NA18563 NA18570 NA18612 NA18621 NA18622 NA18573 NA18577 NA18592 NA18636 NA18593 NA18596 NA18597 NA18536 NA18962 NA18544 NA18546 NA18610 NA18615 NA18557 NA18559 NA18619 NA18638 NA18639 NA18627 NA18628 NA18631 NA18634 NA18745 NA18640 NA18641 NA18748 NA18749 NA18642 NA18599 NA18626 NA18595 NA18618 NA18740
rs11252546 C/T TT TT TT TT TT TT TT TT TT TT TT TT CT TT TT TT CT TT TT TT TT TT TT TT TT TT TT CT TT TT TT TT TT TT TT TT TT TT TT TT CT TT TT TT TT TT TT TT TT TT TT TT TT TT TT TT TT TTTT TT TT TT
rs7909677 A/G AA AA AA AA AG AA AA AA AA AA AA AA AA AA AA AA AA AA AA AA AA AA AA AA AA AG AA AA AA AA AA AA AA AA AG AA AA AA AA AA AA AA AA AA AA AA AA AA AA AA AA AA AA AA AA AG AA AA AG AA AA AA
rs10904494 A/C AA AA AA AA AA AA AA AA AA AA AA AA AC AA AA AA AC AA AA AA AA AA AA AA AA AA AA AC AA AA AA AA AA AA AA AA AA AA AA AA AC AA AA AA AA AA AA AA AA AA AA AA AA AA AA NN AA AAAA NN AA NN
rs11591988 C/T CT TT CC TT CT CT TT CT TT TT TT CT CT CC TT TT CT CT CT TT TT CC TT CC CT CT CT CC TT CC CC CT TT CC CC TT TT CT TT TT CT CT TT CT TT CT TT CT CT TT TT CT CT CT NN TT CC CTCC NN TT NN

$ ll
total 36
drwxr-xr-x  2 hpli Biouser 4096 4月  28 16:44 ./
drwxr-xr-x 10 hpli Biouser 4096 4月  28 10:50 ../
-rw-r--r--  1 hpli Biouser 5116 4月  28 16:42 3-testcase.txt
-rw-r--r--  1 hpli Biouser 5270 4月  28 16:43 3-testcontrol.txt
-rw-r--r--  1 hpli Biouser 9072 4月  28 16:41 4-nucleobase_count.py
$ python2 hapmap_calculation.py 3-testcase.txt 3-testcontrol.txt test.txt test_matrix.txt
SNPID = rs11252546
SNPID = rs7909677
SNPID = rs10904494
SNPID = rs11591988
SNPID = rs4508132
SNPID = rs10904561
SNPID = rs7917054
SNPID = rs7906287
SNPID = rs4495823
SNPID = rs2379076
SNPID = rs11253478
SNPID = rs9419557
SNPID = rs9286070
SNPID = rs9419560
SNPID = rs9419561
SNPID = rs11253562
SNPID = rs4881551
SNPID = rs4880750
SNPID = rs11594819
SNPID = rs9419419
SNPID = rs7909028
SNPID = rs7476951
SNPID = rs12146291
$ ll
total 48
drwxr-xr-x  2 hpli Biouser 4096 4月  28 16:44 ./
drwxr-xr-x 10 hpli Biouser 4096 4月  28 10:50 ../
-rw-r--r--  1 hpli Biouser 5116 4月  28 16:42 3-testcase.txt
-rw-r--r--  1 hpli Biouser 5270 4月  28 16:43 3-testcontrol.txt
-rw-r--r--  1 hpli Biouser 9072 4月  28 16:41 4-nucleobase_count.py
-rw-r--r--  1 hpli Biouser 7440 4月  28 16:44 test_matrix.txt
-rw-r--r--  1 hpli Biouser 1490 4月  28 16:44 test.txt
$ head -5 test_matrix.txt
	Case(1)/Control(0)	rs11252546	rs7909677	rs10904494	rs11591988	rs4508132	rs10904561	rs7917054	rs7906287	rs4495823	rs2379076	rs11253478	rs9419557	rs9286070	rs9419560	rs9419561  rs11253562	rs4881551	rs4880750	rs11594819	rs9419419	rs7909028	rs7476951	rs12146291
NA18532	1	0	0	0	1	1	1	0	1	1	1	1	0	0	0	0	0	1	1	0	0	0	0	1
NA18605	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
NA18542	1	0	0	0	2	2	2	0	2	1	1	2	0	0	0	0	0	1	1	1	0	1	1	1
NA18550	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
$ head -5 test.txt
SNPID	Main_Base	Minor_Base	Case_Main	Case_Minor	Ctr_Main	Ctr_Minor	Tol_#1_AA	Tol_#2_AT	Tol_#3_TT	Tol_NN	Num_of_Sample	Case_#1_AA	Case_#2_AT	Case_#3_TT	Contrl_#1_AA	Contrl_#2_AT	Contrl_#3_TT
rs11252546	T	C	118	4	123	5	116	9	0	1	126	57	4	0	59	5	0
rs7909677	A	G	119	5	121	5	115	10	0	1	126	57	5	0	58	5	0
rs10904494	A	C	112	4	123	5	113	9	0	4	126	54	4	0	59	5	0
rs11591988	T	C	72	46	78	50	50	50	23	3	126	24	24	11	26	26	12  

```
