#!/bin/bash

#This script will generate a pairwise table with all the possible combination between IDs
#After that, it will sum the values of shared IBD segments 
#To speed up the process, it will run in chunks
#Make sure to merge together the different cMmatrices
#To run it in background and output to a log file do the following:
#chmod +x MergeIBD.sh
#./MergeIBD.sh > MergeIBD.log 2>&1 &

thisdir=$(pwd)

#Make sure to unzip all ibd files
gunzip *.ibd.gz

#First, we will paste files from all chromosomes, don't worry about the order
cat *.ibd > all_chr.ibd

#Second, we will generate a file with all possible pariwise comparisons between all IDs
awk '{print $1}' all_chr.ibd | cat > concat.txt
cat all_chr.ibd | awk '{print $3}' >> concat.txt
sort -u concat.txt > concat.sorted.txt

#Generating all possible pairwise combinations without repetition
awk '{
      line[++c] = $1
  }
  END {
      for ( i = 1; i <= c; i++ )
      {
          for ( j = i+1; j <= c; j++ )
          {
              if (line[i] != line[j]) { printf("%s\t%s\n",line[i],line[j]) }
          }
       }
  }
' concat.sorted.txt > combination.txt

#Splitting the combination file, to get manageable number of combinations
split -l 50000 combination.txt combination_split_

counter=0
for file in combination_split_*
do 
	echo "#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=10:00:00
#PBS -l pmem=16gb
#PBS -A open #jlt22_b_g_sc_default or open
#PBS -j oe

#Moving to directory
cd ${thisdir}

awk 'FNR==NR{array[\$1\" \"\$3]+=\$9; barray[\$3\" \"\$1]+=\$9; next} {
	  finarray[\$1\" \"\$2]=0} {for (i in finarray) finarray[i]=barray[i]+array[i]}
	  	 END {for (i in finarray) print i\" \"finarray[i]}' all_chr.ibd ${file} > cMmatrix_${counter}.txt" >> job_${counter}.pbs

  if [ $counter == 99 ]; then
    sleep 30m
  fi

	qsub job_${counter}.pbs
	counter=$((counter+1))
done
