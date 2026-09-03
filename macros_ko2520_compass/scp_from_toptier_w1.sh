#scp top-tier@192.168.1.37:/home/top-tier/ko2520/csi_260901_W1/DAQ/run_370/RAW/* data_compass/
#scp top-tier@192.168.1.37:/home/top-tier/ko2520/csi_260901_W1/DAQ/run_$1/RAW/* data_compass/
rsync -ravzp --progress top-tier@192.168.1.37:/home/top-tier/ko2520/csi_260901_ch16-ch19_CsI/DAQ/run_$1/RAW/* data_compass/
