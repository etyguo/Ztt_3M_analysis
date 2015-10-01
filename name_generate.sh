#!/bin/bash
rm ./python.txt
for i in {1..200}; do
echo "'file:/scratch/osg/etyguo/Ztt_3M_sim/$i.root'," >> ./python.txt
#chmod 744 ./runMyJob.sh
#bsub -q 8nd -J job1 < runMyJob.sh
#head -n -1 ./runMyJob.sh >>temp.sh
#mv temp.sh ./runMyJob.sh 
#echo "$i"
done
