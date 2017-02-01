#!/bin/sh

#first loop
#loopStart,i
for i in `ls -d condition*`; do
    cd $i
    #second loop
    #loopStart,j
    for j in `ls -d sample*`; do
        cd $j
         #third loop
        #loopStart,l
        for l in `ls -f *.txt`; do
            #step1
            #@1,0,copy1: copy original file to .copy1
            cp $l $l.copy1
            #step2
            #@2,1,copy2: copy .copy1 to .copy2
            cp $l.copy1 $l.copy2
        #loopEnd
        done
        #step3
        #@3,1,merge1: merge .copy1 to .copy1.merged
        cat *.copy1 > $j.copy1.merged
        #step4
        #@4,2,merge2: merge .copy2 to .copy2.merged
        cat *.copy2 > $j.copy2.merged
        cd ..
    #loopEnd
    done
    #step5
    #@5,3.4,mergeall: merge everything together
    cat */*merged > $i.all
    cd ..
#loopEnd
done
