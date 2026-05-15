scp ndpsdaq@10.1.100.70:/home/ndpsdaq/online/kobra_daq_2024/output00${1}.root ./data_kobra_daq/
scp ndpsdaq@10.1.100.70:/home/ndpsdaq/online/cens/kobra/root/output00000${1}.root ./data_kobra_ana/
echo lte_tag_${2}.dat
cp ~/lilak/ko2421/macros/raw_data/lte_data/lte_tag_${2}.dat ./data_mte/
