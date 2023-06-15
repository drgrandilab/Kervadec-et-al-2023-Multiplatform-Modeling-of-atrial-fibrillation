run_single_sim() {

	mkdir Bin_$1
	./Atria_2D_MPI_Ghost PLB_Sim $1 Total_time 20000
	mv v_*.bin Bin_$1
}



run_single_sim kmf.1.0.ISO.0.0
run_single_sim kmf.0.25.ISO.0.0
run_single_sim kmf.1.0.ISO.1.0
run_single_sim kmf.0.25.ISO.1.0






# run_APD_BCL () {

# ./GB_Main $1 $2 0 >> APDr.kmf.$2

# mv HAM_wrap_out.dat HAM_wrap_out.dat.BCL.$1.kmf.$2

# }



#  # run_APD_BCL  1 1 1  1 1 1   0.5 0.5 1 365   



# for i in `seq 300 50 1000`
# do 
# echo $i
# run_APD_BCL $i 1.0
# #  run_APD_BCL  0.8054  1.137      0.9191     0.9682     0.9391   1.03435   1 0 1 $i
# done


# for i in `seq 300 50 1000`
# do 
# echo $i
# run_APD_BCL $i 0.5
# #  run_APD_BCL  0.8054  1.137      0.9191     0.9682     0.9391   1.03435   1 0 1 $i
# done




# for i in `seq 300 50 1000`
# do 
# echo $i
# run_APD_BCL $i 0.25
# #  run_APD_BCL  0.8054  1.137      0.9191     0.9682     0.9391   1.03435   1 0 1 $i
# done

