#!/bin/sh

#first loop
for i in `ls -d condition*`; do
    cd $i
    #second loop
    for j in `ls -d sample*`; do
        cd $j
        #third loop
        for l in `ls -f *.txt`; do
            #step1
            cp $l $l.copy1
            #step2
            cp $l.copy1 $l.copy2
        done
        #step3
        cat *.copy1 > $j.copy1.merged
        #step4
        cat *.copy2 > $j.copy2.merged
        cd ..
    done
    #step5
    cat */*merged > $i.all
    cd ..
done
