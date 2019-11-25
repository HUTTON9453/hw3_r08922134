#!/bin/bash

for i in {1..10}
do
	#perl separator_big5.pl test_data/$i.txt > seg_data/$i.txt
	disambig -text seg_data/$i.txt -map ZhuYin-Big5.map -lm bigram.txt -order 2 > result_data/SRILM/bigram/$i.txt
	disambig -text seg_data/$i.txt -map ZhuYin-Big5.map -lm trigram.txt -order 3 > result_data/SRILM/trigram/$i.txt
done
