#!/bin/bash

for i in {1..10}
do
	#perl separator_big5.pl test_data/$i.txt > seg_data/$i.txt
	./mydisambig seg_data/$i.txt ZhuYin-Big5.map bigram.lm result_data/hutton/bigram/$i.txt
done
