#include <stdio.h>
#include "Ngram.h"
#include "VocabMap.h"
#include "LM.h"
#include "File.h"
#include "Vocab.h" 
#include "Prob.h"
#include "iconv.h"



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
void bigram_viterbi(char* output_filename, char* seg_filename, Ngram lm,Vocab voc, VocabMap map, Vocab ZhuYin, Vocab Big5){
	//output File
	FILE *fp;	
	fp = fopen(output_filename,"w+");
	
	char *line;
	File segFile(seg_filename, "r");
	while(line = segFile.getline()){
		//SRILM disambig output <s>X X X</s>
		VocabString sentence[5000];
		unsigned int count = Vocab::parseWords(line, &(sentence[1]), 5000);
		sentence[0] = "<s>";
		sentence[count+1] = "</s>";
		count += 2;
		//build variable for viterbi 
		LogP deltap[200][1024] = {{0.0}};
		VocabIndex VidxGraph[200][1024];
		int Backtrack[200][1024];
		int CandiNum[200];
		Prob p;
		VocabIndex v_idx; 
		VocabIndex empty[] = {Vocab_None};
		VocabIndex bi[] = {Vocab_None, Vocab_None};
		//initialize viterbi
		VocabMapIter iter(map, ZhuYin.getIndex("<s>"));
		iter.init();
		int size = 0;
		while (iter.next(v_idx, p)) {
			VocabIndex wid = voc.getIndex(Big5.getWord(v_idx));
			LogP temp = lm.wordProb(wid, empty);//unigram
			if(temp== LogP_Zero){
				deltap[0][size] = -100;
			}else{
				deltap[0][size] = temp;
			}
			VidxGraph[0][size] = wid;
			Backtrack[0][size] = -1; 
			size++;
		}
		CandiNum[0] = size;
		//recursive solve viterbi
		for (int i = 1; i < count; i++) {
			VocabMapIter iter(map, ZhuYin.getIndex(sentence[i]));
			iter.init();
			size = 0;
			while (iter.next(v_idx, p)) {
				VocabIndex wid = voc.getIndex(Big5.getWord(v_idx));
				if(wid == Vocab_None){
					wid = voc.getIndex(Vocab_Unknown);
				}
				// max(logP(q_i|W_(t-1)))+logW_(1:t-1)
				LogP maxp = LogP_Zero;
				for (int j = 0; j < CandiNum[i-1]; j++) {
						//VocabIndex last = voc.getIndex(Big5.getWord(VidxGraph[i-1][j]));
						VocabIndex last = VidxGraph[i-1][j];
						if(last == Vocab_None){
							last = voc.getIndex(Vocab_Unknown);
						}						
						bi[0] = last;

						LogP logp = lm.wordProb(wid, bi); 
						if (logp == LogP_Zero ){
							logp = -100;
						} 

						logp += deltap[i-1][j];
						if (logp > maxp) {
							maxp = logp;
							Backtrack[i][size] = j;
						}
				}
				deltap[i][size] = maxp;
				VidxGraph[i][size] = wid;
				/*if(deltap[i][size] > 0 && i==count-1){
					finalmaxp = deltap[i][size];
					finalmax = size;
				}*/
				size++;
			}
			CandiNum[i] = size;
		}
		int finalmax = -1,j;
		LogP finalmaxp = LogP_Zero;
		for (j = 0; j < CandiNum[count-1]; j++) {
			if (deltap[count-1][j] > finalmaxp) {
				finalmaxp = deltap[count-1][j];
				finalmax = j;
			}
		}
		VocabString AnsPath[maxWordLength];
		AnsPath[0] = "<s>";
		AnsPath[count-1] = "</s>";
		for (int i = count-1; i > 0; i--) {
			AnsPath[i] = voc.getWord(VidxGraph[i][finalmax]);
			finalmax = Backtrack[i][finalmax];
		}

		// Print the Answer Path
		for (int i = 0; i < count; i++)
			fprintf(fp,"%s%s", AnsPath[i], (i == count-1)? "\n": " ");
	}
	fclose(fp);

}

void trigram_viterbi(char* output_filename, char* seg_filename, Ngram lm,Vocab voc, VocabMap map, Vocab ZhuYin, Vocab Big5){
	//output File
	FILE *fp;	
	fp = fopen(output_filename,"w+");
	
	char *line;
	File segFile(seg_filename, "r");
	while(line = segFile.getline()){
		//SRILM disambig output <s>X X X</s>
		VocabString sentence[500];
		unsigned int count = Vocab::parseWords(line, &(sentence[1]), 500);
		sentence[0] = "<s>";
		sentence[count+1] = "</s>";
		count += 2;
		//build variable for viterbi 
		LogP deltap[200][2000][2000] = {{{0.0}}};
		VocabIndex VidxGraph[200][2000];
		int Backtrack[200][2000][2000];
		int CandiNum[200];
		Prob p;
		VocabIndex v_idx; 
		VocabIndex empty[] = {Vocab_None};
		VocabIndex bi[] = {Vocab_None, Vocab_None};
		VocabIndex tri[] = {Vocab_None, Vocab_None, Vocab_None};
		//initialize viterbi
		VocabMapIter iter(map, ZhuYin.getIndex(sentence[0]));
		iter.init();
		int size = 0;
		while (iter.next(v_idx, p)) {
			VocabIndex wid = voc.getIndex(Big5.getWord(v_idx));
			LogP temp = lm.wordProb(wid, empty);//unigram
			if(temp==LogP_Zero ){
				deltap[0][0][size] = -100;
			}else{
				deltap[0][0][size] = temp;
			}
			VidxGraph[0][size] = v_idx;
			Backtrack[0][0][size] = -1; 
			size++;
		}
		CandiNum[0] = size;
		
		VocabMapIter iter2(map, ZhuYin.getIndex(sentence[1]));
		iter2.init();
		size = 0;
		while (iter2.next(v_idx, p)) {
			VocabIndex wid = voc.getIndex(Big5.getWord(v_idx));
			if(wid == Vocab_None){
				wid = voc.getIndex(Vocab_Unknown);
			}
			// max(logP(q_i|W_(t-1)))+logW_(1:t-1)

			LogP maxp = LogP_Zero;
			for (int j = 0; j < CandiNum[0]; j++) {
				VocabIndex last = voc.getIndex(Big5.getWord(VidxGraph[0][j]));
				if(last == Vocab_None){
					last = voc.getIndex(Vocab_Unknown);
				}						
				bi[0] = last;
				LogP logp = lm.wordProb(wid, bi); 
				if (logp == LogP_Zero ){
					logp = -100;
				} 

				logp += deltap[0][0][0];
				if (logp > maxp) {
					maxp = logp;
					Backtrack[1][0][size] = 0;
				}
			}
			deltap[1][0][size] = maxp;
			VidxGraph[1][size] = v_idx;
			size++;
			
		}
		CandiNum[1] = size;
		//recursive solve viterbi
		for (int i = 2; i < count; i++) {
			VocabMapIter iter(map, ZhuYin.getIndex(sentence[i]));
			iter.init();
			size = 0;
			while (iter.next(v_idx, p)) {
				VocabIndex wid = voc.getIndex(Big5.getWord(v_idx));
				if(wid == Vocab_None){
					wid = voc.getIndex(Vocab_Unknown);
				}
				// max(logP(q_i|q_j,q_k)+logW_(1:t-1)
				for (int j = 0; j < CandiNum[i-1]; j++) {
					VocabIndex last = voc.getIndex(Big5.getWord(VidxGraph[i-1][j]));
					if(last == Vocab_None){
						last = voc.getIndex(Vocab_Unknown);
					}	
					tri[1] = last;					
					LogP maxp = LogP_Zero;
					for(int k = 0;k < CandiNum[i-2];k++){
						VocabIndex last2 = voc.getIndex(Big5.getWord(VidxGraph[i-2][k]));
						if(last2 == Vocab_None){
							last2 = voc.getIndex(Vocab_Unknown);
						}						
						tri[0] = last2;

						LogP logp = lm.wordProb(wid, tri); 
						if (logp == LogP_Zero ){
							logp = -100;
						} 

						logp += deltap[i-1][k][j];
						if (logp > maxp) {
							maxp = logp;
							Backtrack[i][j][size] = k;
						}
					}
					deltap[i][j][size] = maxp;
				}
				VidxGraph[i][size] = v_idx;
				size++;
			}
			CandiNum[i] = size;
		}
		int fmax=0,fmax2=0,qq,j,k=0;
		LogP finalmaxp = -1/0;
		for (j = 0; j < CandiNum[count-2]; j++) {
			if(deltap[count-1][j][0]>finalmaxp){
				finalmaxp = deltap[count-1][j][0];
				fmax = j;
			}
		}
		std::cout<<fmax<<std::endl;
		VocabString AnsPath[100];
		AnsPath[0] = "<s>";
		AnsPath[count-1] = "</s>";
		AnsPath[count-1] = Big5.getWord(VidxGraph[count-1][0]);
		for (int i = count-1; i > 0; i--) {
			//AnsPath[i] = Big5.getWord(VidxGraph[i][fmax2]);
			qq = Backtrack[i][fmax][fmax2];
			fmax = fmax2;
			fmax2 = qq;
		}
		// Print the Answer Path
	//	for (int i = 0; i < count; i++)
	//		fprintf(fp,"%s%s", AnsPath[i], (i == count-1)? "\n": " ");
	}
	fclose(fp);

}

int main(int argc, char *argv[]){
	
	int ngram_order;
	char* seg_filename;
	char* map_filename; 
	char* lm_filename;
	char* output_filename;
	if(argc==5){
		ngram_order = 2;
		seg_filename = argv[1];
		map_filename = argv[2]; 
		lm_filename = argv[3];
		output_filename =argv[4];
	}else if(argc==6){
		ngram_order = atoi(argv[5]);
		seg_filename = argv[1];
		map_filename = argv[2]; 
		lm_filename = argv[3];
		output_filename =argv[4];
	}else{
		printf("two types commands: \n");	
		printf("./mydisambig $seg_filename $map_filename $LM_filename $output_filename\n");
		printf("./mydisambig $seg_filename $map_filename $LM_filename $output_filename $ngram_order\n");
		exit(1);
	}
	
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

	if(ngram_order==2){
		bigram_viterbi(output_filename, seg_filename, lm, voc, map, ZhuYin, Big5);
	}else{
		trigram_viterbi(output_filename, seg_filename, lm, voc, map, ZhuYin, Big5);
	}

	return 0;
}
