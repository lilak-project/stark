##################################
## e4 ############################
##################################
#for id in $(seq 20 38)
#do
#    RUN=$id root run_reco.C
#done
#
###################################
### e5 ############################
###################################
#for id in $(seq 39 57)
#do
#    RUN=$id root run_reco.C
#done
#
#for id in $(seq 59 72)
#do
#    RUN=$id root run_reco.C
#done
#
##################################
## e8 ############################
##################################
for id in $(seq 81 83)
do
    RUN=$id root run_reco.C
done
#RUN=92 root run_reco.C
#for id in $(seq 94 97)
#do
#    RUN=$id root run_reco.C
#done
for id in $(seq 99 107)
do
    RUN=$id root run_reco.C
done
for id in $(seq 110 111)
do
    RUN=$id root run_reco.C
done
##################################
