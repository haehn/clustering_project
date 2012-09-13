# bsub -o /dev/null -e error.%J.%I -J "nlv1[1-50]" bash ../clustering_project/tempdir_wrapper.sh python doclustering.py -d ~/storage/regime2/nni/level3/sim
# sleep 0.5
# bsub -o /dev/null -e error.%J.%I -J "nlv1[1-50]" bash ../clustering_project/tempdir_wrapper.sh python doclustering.py -d ~/storage/regime2/nni/level4/sim
# sleep 0.5
# bsub -o /dev/null -e error.%J.%I -J "nlv1[1-50]" bash ../clustering_project/tempdir_wrapper.sh python doclustering.py -d ~/storage/regime2/nni/level5/sim
# sleep 0.5
# bsub -o /dev/null -e error.%J.%I -J "nlv1[1-50]" bash ../clustering_project/tempdir_wrapper.sh python doclustering.py -d ~/storage/regime2/nni/level6/sim
# sleep 0.5
# bsub -o /dev/null -e error.%J.%I -J "nlv1[1-50]" bash ../clustering_project/tempdir_wrapper.sh python doclustering.py -d ~/storage/regime2/nni/level7/sim
# sleep 0.5
# bsub -o /dev/null -e error.%J.%I -J "nlv1[1-50]" bash ../clustering_project/tempdir_wrapper.sh python doclustering.py -d ~/storage/regime2/nni/level8/sim
# sleep 0.5
# bsub -o /dev/null -e error.%J.%I -J "nlv1[1-50]" bash ../clustering_project/tempdir_wrapper.sh python doclustering.py -d ~/storage/regime2/nni/level9/sim
# sleep 0.5
# bsub -o /dev/null -e error.%J.%I -J "nlv1[1-50]" bash ../clustering_project/tempdir_wrapper.sh python doclustering.py -d ~/storage/regime2/nni/level10/sim
# sleep 0.5
# bsub -o /dev/null -e error.%J.%I -J "slv1[1-50]" bash ../clustering_project/tempdir_wrapper.sh python doclustering.py -d ~/storage/regime2/spr/level1/sim
# sleep 0.5
# bsub -o /dev/null -e error.%J.%I -J "slv1[1-50]" bash ../clustering_project/tempdir_wrapper.sh python doclustering.py -d ~/storage/regime2/spr/level2/sim 
# sleep 0.5
# bsub -o /dev/null -e error.%J.%I -J "slv1[1-50]" bash ../clustering_project/tempdir_wrapper.sh python doclustering.py -d ~/storage/regime2/spr/level3/sim
# sleep 0.5
# bsub -o /dev/null -e error.%J.%I -J "slv1[1-50]" bash ../clustering_project/tempdir_wrapper.sh python doclustering.py -d ~/storage/regime2/spr/level4/sim
# sleep 0.5
# bsub -o /dev/null -e error.%J.%I -J "slv1[1-50]" bash ../clustering_project/tempdir_wrapper.sh python doclustering.py -d ~/storage/regime2/spr/level5/sim
# sleep 0.5
# bsub -o /dev/null -e error.%J.%I -J "slv1[1-50]" bash ../clustering_project/tempdir_wrapper.sh python doclustering.py -d ~/storage/regime2/spr/level6/sim
# sleep 0.5
# bsub -o /dev/null -e error.%J.%I -J "slv1[1-50]" bash ../clustering_project/tempdir_wrapper.sh python doclustering.py -d ~/storage/regime2/spr/level7/sim
# sleep 0.5
# bsub -o /dev/null -e error.%J.%I -J "slv1[1-50]" bash ../clustering_project/tempdir_wrapper.sh python doclustering.py -d ~/storage/regime2/spr/level8/sim
# sleep 0.5
# bsub -o /dev/null -e error.%J.%I -J "slv1[1-50]" bash ../clustering_project/tempdir_wrapper.sh python doclustering.py -d ~/storage/regime2/spr/level9/sim
# sleep 0.5
# bsub -o /dev/null -e error.%J.%I -J "slv1[1-50]" bash ../clustering_project/tempdir_wrapper.sh python doclustering.py -d ~/storage/regime2/spr/level10/sim
# sleep 0.5
# bsub -o /dev/null -e error.%J.%I -J "clv1[1-50]" bash ../clustering_project/tempdir_wrapper.sh python doclustering.py -d ~/storage/regime2/coal/level1/sim
# sleep 0.5
# bsub -o /dev/null -e error.%J.%I -J "clv1[1-50]" bash ../clustering_project/tempdir_wrapper.sh python doclustering.py -d ~/storage/regime2/coal/level2/sim
# sleep 0.5
# bsub -o /dev/null -e error.%J.%I -J "clv1[1-50]" bash ../clustering_project/tempdir_wrapper.sh python doclustering.py -d ~/storage/regime2/coal/level4/sim
# sleep 0.5
# bsub -o /dev/null -e error.%J.%I -J "clv1[1-50]" bash ../clustering_project/tempdir_wrapper.sh python doclustering.py -d ~/storage/regime2/coal/level8/sim
# sleep 0.5
# bsub -o /dev/null -e error.%J.%I -J "clv1[1-50]" bash ../clustering_project/tempdir_wrapper.sh python doclustering.py -d ~/storage/regime2/coal/level16/sim
# sleep 0.5
# bsub -o /dev/null -e error.%J.%I -J "clv1[1-50]" bash ../clustering_project/tempdir_wrapper.sh python doclustering.py -d ~/storage/regime2/coal/level32/sim
# sleep 0.5
# bsub -o /dev/null -e error.%J.%I -J "clv1[1-50]" bash ../clustering_project/tempdir_wrapper.sh python doclustering.py -d ~/storage/regime2/coal/level64/sim
# sleep 0.5
# bsub -o /dev/null -e error.%J.%I -J "clv1[1-50]" bash ../clustering_project/tempdir_wrapper.sh python doclustering.py -d ~/storage/regime2/coal/level128/sim
# sleep 0.5



# BATCH
# BIONJ
# COMMANDS

# bsub -o /dev/null -J "njn1[1-50]" bash ../clustering_project/tempdir_wrapper.sh python runbatchbionj.py -d ~/storage/regime2/nni/level1/sim -p bionj_clustering
# sleep 0.5
# bsub -o /dev/null -J "njn2[1-50]" bash ../clustering_project/tempdir_wrapper.sh python runbatchbionj.py -d ~/storage/regime2/nni/level2/sim -p bionj_clustering
# sleep 0.5
# bsub -o /dev/null -J "njn3[1-50]" bash ../clustering_project/tempdir_wrapper.sh python runbatchbionj.py -d ~/storage/regime2/nni/level3/sim -p bionj_clustering
# sleep 0.5
# bsub -o /dev/null -J "njn4[1-50]" bash ../clustering_project/tempdir_wrapper.sh python runbatchbionj.py -d ~/storage/regime2/nni/level4/sim -p bionj_clustering
# sleep 0.5
# bsub -o /dev/null -J "njn5[1-50]" bash ../clustering_project/tempdir_wrapper.sh python runbatchbionj.py -d ~/storage/regime2/nni/level5/sim -p bionj_clustering
# sleep 0.5
# bsub -o /dev/null -J "njn6[1-50]" bash ../clustering_project/tempdir_wrapper.sh python runbatchbionj.py -d ~/storage/regime2/nni/level6/sim -p bionj_clustering
# sleep 0.5
# bsub -o /dev/null -J "njn7[1-50]" bash ../clustering_project/tempdir_wrapper.sh python runbatchbionj.py -d ~/storage/regime2/nni/level7/sim -p bionj_clustering
# sleep 0.5
# bsub -o /dev/null -J "njn8[1-50]" bash ../clustering_project/tempdir_wrapper.sh python runbatchbionj.py -d ~/storage/regime2/nni/level8/sim -p bionj_clustering
# sleep 0.5
# bsub -o /dev/null -J "njn9[1-50]" bash ../clustering_project/tempdir_wrapper.sh python runbatchbionj.py -d ~/storage/regime2/nni/level9/sim -p bionj_clustering
# sleep 0.5
# bsub -o /dev/null -J "njn10[1-50]" bash ../clustering_project/tempdir_wrapper.sh python runbatchbionj.py -d ~/storage/regime2/nni/level10/sim -p bionj_clustering
# sleep 0.5
# bsub -o /dev/null -J "njs1[1-50]" bash ../clustering_project/tempdir_wrapper.sh python runbatchbionj.py -d ~/storage/regime2/spr/level1/sim -p bionj_clustering
# sleep 0.5
# bsub -o /dev/null -J "njs2[1-50]" bash ../clustering_project/tempdir_wrapper.sh python runbatchbionj.py -d ~/storage/regime2/spr/level2/sim -p bionj_clustering
# sleep 0.5
# bsub -o /dev/null -J "njs3[1-50]" bash ../clustering_project/tempdir_wrapper.sh python runbatchbionj.py -d ~/storage/regime2/spr/level3/sim -p bionj_clustering
# sleep 0.5
# bsub -o /dev/null -J "njs4[1-50]" bash ../clustering_project/tempdir_wrapper.sh python runbatchbionj.py -d ~/storage/regime2/spr/level4/sim -p bionj_clustering
# sleep 0.5
# bsub -o /dev/null -J "njs5[1-50]" bash ../clustering_project/tempdir_wrapper.sh python runbatchbionj.py -d ~/storage/regime2/spr/level5/sim -p bionj_clustering
# sleep 0.5
# bsub -o /dev/null -J "njs6[1-50]" bash ../clustering_project/tempdir_wrapper.sh python runbatchbionj.py -d ~/storage/regime2/spr/level6/sim -p bionj_clustering
# sleep 0.5
# bsub -o /dev/null -J "njs7[1-50]" bash ../clustering_project/tempdir_wrapper.sh python runbatchbionj.py -d ~/storage/regime2/spr/level7/sim -p bionj_clustering
# sleep 0.5
# bsub -o /dev/null -J "njs8[1-50]" bash ../clustering_project/tempdir_wrapper.sh python runbatchbionj.py -d ~/storage/regime2/spr/level8/sim -p bionj_clustering
# sleep 0.5
# bsub -o /dev/null -J "njs9[1-50]" bash ../clustering_project/tempdir_wrapper.sh python runbatchbionj.py -d ~/storage/regime2/spr/level9/sim -p bionj_clustering
# sleep 0.5
# bsub -o /dev/null -J "njs10[1-50]" bash ../clustering_project/tempdir_wrapper.sh python runbatchbionj.py -d ~/storage/regime2/spr/level10/sim -p bionj_clustering
# sleep 0.5
# bsub -o /dev/null -J "njc1[1-50]" bash ../clustering_project/tempdir_wrapper.sh python runbatchbionj.py -d ~/storage/regime2/coal/level1/sim -p bionj_clustering
# sleep 0.5
# bsub -o /dev/null -J "njc2[1-50]" bash ../clustering_project/tempdir_wrapper.sh python runbatchbionj.py -d ~/storage/regime2/coal/level2/sim -p bionj_clustering
# sleep 0.5
# bsub -o /dev/null -J "njc3[1-50]" bash ../clustering_project/tempdir_wrapper.sh python runbatchbionj.py -d ~/storage/regime2/coal/level4/sim -p bionj_clustering
# sleep 0.5
# bsub -o /dev/null -J "njc4[1-50]" bash ../clustering_project/tempdir_wrapper.sh python runbatchbionj.py -d ~/storage/regime2/coal/level8/sim -p bionj_clustering
# sleep 0.5
# bsub -o /dev/null -J "njc5[1-50]" bash ../clustering_project/tempdir_wrapper.sh python runbatchbionj.py -d ~/storage/regime2/coal/level16/sim -p bionj_clustering
# sleep 0.5
# bsub -o /dev/null -J "njc6[1-50]" bash ../clustering_project/tempdir_wrapper.sh python runbatchbionj.py -d ~/storage/regime2/coal/level32/sim -p bionj_clustering
# sleep 0.5
# bsub -o /dev/null -J "njc7[1-50]" bash ../clustering_project/tempdir_wrapper.sh python runbatchbionj.py -d ~/storage/regime2/coal/level64/sim -p bionj_clustering
# sleep 0.5
# bsub -o /dev/null -J "njc8[1-50]" bash ../clustering_project/tempdir_wrapper.sh python runbatchbionj.py -d ~/storage/regime2/coal/level128/sim -p bionj_clustering
# sleep 0.5

# BATCH
# PHYML
# COMMANDS

# bsub -o /dev/null -J "slv1[1-1000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level1/1-1000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "slv1[1001-2000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level1/1001-2000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "slv1[2001-3000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level1/2001-3000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "slv1[3001-4000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level1/3001-4000/phylip.
# sleep 0.5
bsub -o /dev/null -J "slv1[4001-4400]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level1/4001-5000/phylip.
sleep 0.5

# bsub -o /dev/null -J "slv2[1-1000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level2/1-1000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "slv2[1001-2000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level2/1001-2000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "slv2[2001-3000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level2/2001-3000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "slv2[3001-4000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level2/3001-4000/phylip.
# sleep 0.5
bsub -o /dev/null -J "slv2[4001-4400]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level2/4001-5000/phylip.
sleep 0.5

# bsub -o /dev/null -J "slv3[1-1000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level3/1-1000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "slv3[1001-2000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level3/1001-2000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "slv3[2001-3000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level3/2001-3000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "slv3[3001-4000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level3/3001-4000/phylip.
# sleep 0.5
bsub -o /dev/null -J "slv3[4001-4400]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level3/4001-5000/phylip.
sleep 0.5

# bsub -o /dev/null -J "slv4[1-1000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level4/1-1000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "slv4[1001-2000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level4/1001-2000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "slv4[2001-3000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level4/2001-3000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "slv4[3001-4000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level4/3001-4000/phylip.
# sleep 0.5
bsub -o /dev/null -J "slv4[4001-4400]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level4/4001-5000/phylip.
sleep 0.5

# bsub -o /dev/null -J "slv5[1-1000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level5/1-1000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "slv5[1001-2000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level5/1001-2000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "slv5[2001-3000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level5/2001-3000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "slv5[3001-4000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level5/3001-4000/phylip.
# sleep 0.5
bsub -o /dev/null -J "slv5[4001-4400]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level5/4001-5000/phylip.
sleep 0.5

# bsub -o /dev/null -J "slv6[1-1000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level6/1-1000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "slv6[1001-2000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level6/1001-2000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "slv6[2001-3000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level6/2001-3000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "slv6[3001-4000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level6/3001-4000/phylip.
# sleep 0.5
bsub -o /dev/null -J "slv6[4001-4400]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level6/4001-5000/phylip.
sleep 0.5

# bsub -o /dev/null -J "slv7[1-1000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level7/1-1000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "slv7[1001-2000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level7/1001-2000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "slv7[2001-3000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level7/2001-3000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "slv7[3001-4000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level7/3001-4000/phylip.
# sleep 0.5
bsub -o /dev/null -J "slv7[4001-4400]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level7/4001-5000/phylip.
sleep 0.5

# bsub -o /dev/null -J "slv8[1-1000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level8/1-1000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "slv8[1001-2000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level8/1001-2000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "slv8[2001-3000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level8/2001-3000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "slv8[3001-4000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level8/3001-4000/phylip.
# sleep 0.5
bsub -o /dev/null -J "slv8[4001-4400]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level8/4001-5000/phylip.
sleep 0.5

# bsub -o /dev/null -J "slv9[1-1000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level9/1-1000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "slv9[1001-2000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level9/1001-2000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "slv9[2001-3000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level9/2001-3000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "slv9[3001-4000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level9/3001-4000/phylip.
# sleep 0.5
bsub -o /dev/null -J "slv9[4001-4400]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level9/4001-5000/phylip.
sleep 0.5

# bsub -o /dev/null -J "slv10[1-1000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level10/1-1000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "slv10[1001-2000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level10/1001-2000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "slv10[2001-3000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level10/2001-3000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "slv10[3001-4000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level10/3001-4000/phylip.
# sleep 0.5
bsub -o /dev/null -J "slv10[4001-4400]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/spr/level10/4001-5000/phylip.
sleep 0.5

# bsub -o /dev/null -J "nlv1[1-1000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level1/1-1000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "nlv1[1001-2000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level1/1001-2000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "nlv1[2001-3000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level1/2001-3000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "nlv1[3001-4000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level1/3001-4000/phylip.
# sleep 0.5
bsub -o /dev/null -J "nlv1[4001-4400]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level1/4001-5000/phylip.
sleep 0.5

# bsub -o /dev/null -J "nlv2[1-1000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level2/1-1000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "nlv2[1001-2000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level2/1001-2000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "nlv2[2001-3000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level2/2001-3000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "nlv2[3001-4000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level2/3001-4000/phylip.
# sleep 0.5
bsub -o /dev/null -J "nlv2[4001-4400]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level2/4001-5000/phylip.
sleep 0.5

# bsub -o /dev/null -J "nlv3[1-1000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level3/1-1000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "nlv3[1001-2000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level3/1001-2000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "nlv3[2001-3000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level3/2001-3000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "nlv3[3001-4000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level3/3001-4000/phylip.
# sleep 0.5
bsub -o /dev/null -J "nlv3[4001-4400]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level3/4001-5000/phylip.
sleep 0.5

# bsub -o /dev/null -J "nlv4[1-1000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level4/1-1000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "nlv4[1001-2000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level4/1001-2000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "nlv4[2001-3000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level4/2001-3000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "nlv4[3001-4000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level4/3001-4000/phylip.
# sleep 0.5
bsub -o /dev/null -J "nlv4[4001-4400]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level4/4001-5000/phylip.
sleep 0.5

# bsub -o /dev/null -J "nlv5[1-1000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level5/1-1000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "nlv5[1001-2000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level5/1001-2000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "nlv5[2001-3000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level5/2001-3000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "nlv5[3001-4000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level5/3001-4000/phylip.
# sleep 0.5
bsub -o /dev/null -J "nlv5[4001-4400]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level5/4001-5000/phylip.
sleep 0.5

# bsub -o /dev/null -J "nlv6[1-1000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level6/1-1000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "nlv6[1001-2000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level6/1001-2000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "nlv6[2001-3000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level6/2001-3000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "nlv6[3001-4000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level6/3001-4000/phylip.
# sleep 0.5
bsub -o /dev/null -J "nlv6[4001-4400]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level6/4001-5000/phylip.
sleep 0.5

# bsub -o /dev/null -J "nlv7[1-1000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level7/1-1000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "nlv7[1001-2000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level7/1001-2000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "nlv7[2001-3000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level7/2001-3000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "nlv7[3001-4000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level7/3001-4000/phylip.
# sleep 0.5
bsub -o /dev/null -J "nlv7[4001-4400]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level7/4001-5000/phylip.
sleep 0.5

# bsub -o /dev/null -J "nlv8[1-1000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level8/1-1000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "nlv8[1001-2000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level8/1001-2000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "nlv8[2001-3000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level8/2001-3000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "nlv8[3001-4000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level8/3001-4000/phylip.
# sleep 0.5
bsub -o /dev/null -J "nlv8[4001-4400]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level8/4001-5000/phylip.
sleep 0.5

# bsub -o /dev/null -J "nlv9[1-1000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level9/1-1000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "nlv9[1001-2000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level9/1001-2000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "nlv9[2001-3000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level9/2001-3000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "nlv9[3001-4000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level9/3001-4000/phylip.
# sleep 0.5
bsub -o /dev/null -J "nlv9[4001-4400]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level9/4001-5000/phylip.
sleep 0.5

# bsub -o /dev/null -J "nlv10[1-1000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level10/1-1000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "nlv10[1001-2000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level10/1001-2000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "nlv10[2001-3000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level10/2001-3000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "nlv10[3001-4000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level10/3001-4000/phylip.
# sleep 0.5
bsub -o /dev/null -J "nlv10[4001-4400]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/nni/level10/4001-5000/phylip.
sleep 0.5

# bsub -o /dev/null -J "clv1[1-1000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/coal/level1/1-1000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "clv1[1001-2000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/coal/level1/1001-2000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "clv1[2001-3000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/coal/level1/2001-3000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "clv1[3001-4000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/coal/level1/3001-4000/phylip.
# sleep 0.5
bsub -o /dev/null -J "clv1[4001-4400]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/coal/level1/4001-5000/phylip.
sleep 0.5

# bsub -o /dev/null -J "clv2[1-1000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/coal/level2/1-1000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "clv2[1001-2000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/coal/level2/1001-2000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "clv2[2001-3000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/coal/level2/2001-3000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "clv2[3001-4000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/coal/level2/3001-4000/phylip.
# sleep 0.5
bsub -o /dev/null -J "clv2[4001-4400]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/coal/level2/4001-5000/phylip.
sleep 0.5

# bsub -o /dev/null -J "clv4[1-1000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/coal/level4/1-1000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "clv4[1001-2000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/coal/level4/1001-2000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "clv4[2001-3000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/coal/level4/2001-3000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "clv4[3001-4000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/coal/level4/3001-4000/phylip.
# sleep 0.5
bsub -o /dev/null -J "clv4[4001-4400]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/coal/level4/4001-5000/phylip.
sleep 0.5

# bsub -o /dev/null -J "clv8[1-1000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/coal/level8/1-1000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "clv8[1001-2000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/coal/level8/1001-2000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "clv8[2001-3000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/coal/level8/2001-3000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "clv8[3001-4000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/coal/level8/3001-4000/phylip.
# sleep 0.5
bsub -o /dev/null -J "clv8[4001-4400]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/coal/level8/4001-5000/phylip.
sleep 0.5

# bsub -o /dev/null -J "clv16[1-1000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/coal/level16/1-1000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "clv16[1001-2000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/coal/level16/1001-2000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "clv16[2001-3000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/coal/level16/2001-3000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "clv16[3001-4000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/coal/level16/3001-4000/phylip.
# sleep 0.5
bsub -o /dev/null -J "clv16[4001-4400]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/coal/level16/4001-5000/phylip.
sleep 0.5

# bsub -o /dev/null -J "clv32[1-1000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/coal/level32/1-1000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "clv32[1001-2000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/coal/level32/1001-2000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "clv32[2001-3000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/coal/level32/2001-3000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "clv32[3001-4000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/coal/level32/3001-4000/phylip.
# sleep 0.5
bsub -o /dev/null -J "clv32[4001-4400]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/coal/level32/4001-5000/phylip.
sleep 0.5

# bsub -o /dev/null -J "clv64[1-1000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/coal/level64/1-1000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "clv64[1001-2000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/coal/level64/1001-2000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "clv64[2001-3000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/coal/level64/2001-3000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "clv64[3001-4000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/coal/level64/3001-4000/phylip.
# sleep 0.5
bsub -o /dev/null -J "clv64[4001-4400]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/coal/level64/4001-5000/phylip.
sleep 0.5

# bsub -o /dev/null -J "clv128[1-1000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/coal/level128/1-1000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "clv128[1001-2000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/coal/level128/1001-2000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "clv128[2001-3000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/coal/level128/2001-3000/phylip.
# sleep 0.5
# bsub -o /dev/null -J "clv128[3001-4000]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/coal/level128/3001-4000/phylip.
# sleep 0.5
bsub -o /dev/null -J "clv128[4001-4400]" bash ../clustering_project/tempdir_wrapper.sh python runbatchphyml.py -f ~/storage/regime2/coal/level128/4001-5000/phylip.
sleep 0.5
