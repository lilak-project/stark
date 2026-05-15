for id in 199 168 303
do
    RUN=$id root -q -b -l make_c3_histograms.C\(0\)
done

#for id in 199 168 303
#do
#    RUN=$id root -q -b -l make_histograms.C\(0\)
#    RUN=$id root -q -b -l run_slope_correction.C\(0\)
#    RUN=$id root -q -b -l run_energy_calibration.C\(0\)
#    RUN=$id root -q -b -l make_c1_histograms.C\(0\)
#    RUN=$id root -q -b -l run_ballistic_correction.C\(0\)
#    RUN=$id root -q -b -l make_c2_histograms.C\(0\)
#    RUN=$id root -q -b -l run_junction_energy_calibration.C\(0\)
#    RUN=$id root -q -b -l make_c3_histograms.C\(0\)
#done
#
#root -q -b -l run_recalibration_168_with_199_parameters.C
#root -q -b -l summary.C\(0\)

#for id in 303 113 199 168
#do
#    #RUN=$id root -q -b -l make_histograms.C\(0\)
#    #RUN=$id root -q -b -l run_slope_correction.C\(0\)
#    #RUN=$id root -q -b -l run_energy_calibration.C\(0\)
#    #RUN=$id root -q -b -l make_c1_histograms.C\(0\)
#    #RUN=$id root -q -b -l run_ballistic_correction.C\(0\)
#    #RUN=$id root -q -b -l make_c2_histograms.C\(0\)
#    #RUN=$id root -q -b -l run_position_correction.C\(0\)
#    RUN=$id root -q -b -l make_c3_histograms.C\(0\)
#done
