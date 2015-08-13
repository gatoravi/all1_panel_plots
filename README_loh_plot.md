#The JIRA for this issue is here - https://jira.gsc.wustl.edu/browse/TD-929
#The model-group for the 4 clin-seq models(relapse1,2 primary 1,2)is 3e86c50b5b314612a3a76cdfebbccec4

cp /gscmnt/gc13029/info/model_data/950083c2e3bd431b8ca65cbf7b202c9c/buildd6a324a1034f4729a87cd16f298a42bb/AML103/loh/loh.infile primary1_loh.infile
cp /gscmnt/gc13029/info/model_data/950083c2e3bd431b8ca65cbf7b202c9c/buildd6a324a1034f4729a87cd16f298a42bb/AML103/loh/loh.segments.cbs primary1_segments.cbs

cp /gscmnt/gc13029/info/model_data/6ee6807bc01944b78ff2475935dc1903/build7c0532e4606f4f6a85721ebf289d0b3d/AML103/loh/loh.infile primary2_loh.infile
cp /gscmnt/gc13029/info/model_data/6ee6807bc01944b78ff2475935dc1903/build7c0532e4606f4f6a85721ebf289d0b3d/AML103/loh/loh.segments.cbs primary2_segments.cbs

cp /gscmnt/gc13027/info/model_data/c82164783bab4bc98c9402609a4a8221/build6b0915e3dea8425086f46cbd819f5974/AML103/loh/loh.infile relapse1_loh.infile
cp /gscmnt/gc13027/info/model_data/c82164783bab4bc98c9402609a4a8221/build6b0915e3dea8425086f46cbd819f5974/AML103/loh/loh.segments.cbs relapse1_segments.cbs

cp /gscmnt/gc13027/info/model_data/52f3b8ad88fa4df79196d179aa29e00b/build4c73477e164b470eab1348c90396b435/AML103/loh/loh.infile relapse2_loh.infile
cp /gscmnt/gc13027/info/model_data/52f3b8ad88fa4df79196d179aa29e00b/build4c73477e164b470eab1348c90396b435/AML103/loh/loh.segments.cbs relapse2_segments.cbs

#Join all the infiles together.
rm -f all_loh.infile; python  ../../src/join_lohinfile.py primary1_loh.infile primary2_loh.infile relapse1_loh.infile relapse2_loh.infile | sort -k1,1 -nk2,2 | grep -v chrom > all_loh.infile
