#include <iostream>
#include <stdio.h>
#include "Ngram.h"
#include "VocabMap.h"
#include "LM.h"
#include "File.h"
#include "Vocab.h" 
#include "Prob.h"
#include "iconv.h"
using namespace std;
#ifndef CANDIDATE
#define CANDIDATE 1000
#endif

#ifndef LEN
#define LEN 200
#endif

#ifndef SMALL
#define SMALL -100
#endif

void big5toutf8(char* inbuf){
	iconv_t cd = iconv_open("UTF-8","BIG5");	
	char outbuf[1024];
	char *tmpin = inbuf;
  	char *tmpout = outbuf;
  	size_t insize = strlen(tmpin);
   	size_t outsize = 1024;
	size_t ret = iconv (cd, &tmpin, &insize, &tmpout, &outsize);
	printf("inbuf=%s, outbuf=%s\n", inbuf, outbuf);
}

int main(int argc, char *argv[]){
	int ngram_order = 2;
	const char* seg_filename = argv[1];
	const char* map_filename = argv[2]; 
	const char* lm_filename = argv[3];
	const char* output_filename =argv[4];
	
	/*output File*/
	FILE *fp;	
	fp = fopen(output_filename,"w+");
	/*load LM*/
	Vocab voc;
	Ngram lm(voc, ngram_order);
	File lmFile( lm_filename, "r" );
        lm.read(lmFile);
        lmFile.close();
	/*load Map*/
	Vocab ZhuYin, Big5;
	VocabMap map(ZhuYin, Big5);
	File mapFile(map_filename, "r");
	map.read(mapFile);
	mapFile.close();

	/*run all text*/
  	char *line;
	File textFile("./seg_data/1.txt", "r");
  	while(line = textFile.getline()){
		/*SRILM disambig output <s>X X X</s>*/
		VocabString sen[maxWordsPerLine];
		unsigned int count = Vocab::parseWords(line, &(sen[1]), maxWordsPerLine);
		sen[0] = "<s>";
		sen[count+1] = "</s>";
		count += 2;
		/*build variable for viterbi */
		LogP Proba[LEN][1024] = {{0.0}};
		VocabIndex VidxGraph[LEN][1024];
		int Backtrack[LEN][1024];
		int CandiNum[LEN];

		Prob p;
		VocabIndex v_idx; 
		VocabIndex empty_context[] = {Vocab_None};
		VocabIndex bi_context[] = {Vocab_None, Vocab_None};
		/*initialize viterbi*/
		VocabMapIter iter(map, ZhuYin.getIndex(sen[0]));
		iter.init();
		int size = 0;
		while (iter.next(v_idx, p)) {
			VocabIndex candi = voc.getIndex(Big5.getWord(v_idx));
			LogP logp = lm.wordProb(candi, empty_context); //unigram for start <s>
			Proba[0][size] = (logp == -1.0/0.0)? SMALL: logp;
			//printf("%f %f\n",logp, LogP_Zero);
			VidxGraph[0][size] = v_idx;
			Backtrack[0][size] = -1; // start: -1
			size++;
		}
		CandiNum[0] = size;
		/*recursive solve viterbi*/
		for (int i = 1; i < count; i++) {
			VocabMapIter iter(map, ZhuYin.getIndex(sen[i]));
			iter.init();
			size = 0;
			while (iter.next(v_idx, p)) {
				VocabIndex candi = voc.getIndex(Big5.getWord(v_idx));
				candi = (candi == Vocab_None)? voc.getIndex(Vocab_Unknown): candi;

				// See last column which has the highest probability
				LogP maxp = LogP_Zero;
				for (int j = 0; j < CandiNum[i-1]; j++) {
						VocabIndex last = voc.getIndex(Big5.getWord(VidxGraph[i-1][j]));
						last = (last == Vocab_None)? voc.getIndex(Vocab_Unknown): last;
						bi_context[0] = last;

						LogP logp = lm.wordProb(candi, bi_context);
						// Check backoff!!! VERY IMPORTANT!!! 
						//LogP backoff = lm.wordProb(candi, empty_context);
						//if (logp == LogP_Zero && backoff == LogP_Zero) logp = SMALL;
						if (logp == LogP_Zero ) logp = SMALL;

						logp += Proba[i-1][j];
						if (logp > maxp) {
							maxp = logp;
							Backtrack[i][size] = j;
						}
				}
				Proba[i][size] = maxp;
				VidxGraph[i][size] = v_idx;
				size++;
			}
			CandiNum[i] = size;
		}
		LogP maxp = LogP_Zero;
		int max_col = -1, j;
		for (j = 0; j < CandiNum[count-1]; j++) {
			if (Proba[count-1][j] > maxp) {
				maxp = Proba[count-1][j];
				max_col = j;
			}
		}

		printf("%d\n",max_col);
		VocabString AnsPath[maxWordLength];
		AnsPath[0] = "<s>";
		AnsPath[count-1] = "</s>";
		for (int i = count-1; i > 0; i--) {
			AnsPath[i] = Big5.getWord(VidxGraph[i][max_col]);
			max_col = Backtrack[i][max_col];
		}

		// Print the Answer Path
		for (int i = 0; i < count; i++)
			fprintf(fp,"%s%s", AnsPath[i], (i == count-1)? "\n": " ");
	}
	fclose(fp);


	return 0;
}
