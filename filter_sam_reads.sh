#!/bin/bash

/usr/share/java-1.8.0/jdk1.8.0_92/bin/java -Xmx4g -jar /usr/share/picard/picard-tools-2.5.0/picard.jar FilterSamReads I=$1 O='NO_MT3243_v501_1169_FEMALE.sam' READ_LIST_FILE=$2 FILTER=excludeReadList
