#Merge all the cnvs.hq
rm -f all_cnvs.hq; python join_cnvshq.py primary1_cnvs.hq  primary2_cnvs.hq relapse1_cnvs.hq relapse2_cnvs.hq | sort -k1,1 -nk2,2 | sed '1d' > all_cnvs.hq

#Plot the CNVs
rm -f ~/test_p2.pdf; Rscript-3.1.3 plot_cnv_panel.R all_cnvs.hq
