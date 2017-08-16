#!/bin/bash

/usr/share/java-1.8.0/jdk1.8.0_92/bin/java -Xmx4g -jar /usr/share/picard/picard-tools-2.5.0/picard.jar SortSam I=$1 O='sorted.bam' SORT_ORDER=coordinate
