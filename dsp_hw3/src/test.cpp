#include <stdio.h>
#include "Ngram.h"
#include "File.h"
#include "Vocab.h" 
#include "Prob.h"
#include "iconv.h"
#include <map>
#include <vector>
#include <string>

struct cmp_str
{
    bool operator()(char const *a, char const *b)
    {
        return strcmp(a,b) < 0;
    }
};

int parseWords(char *sentence,char **words, int max)
{
  char *word;
  int i = 0;

  const char *const wordSeparators = " \t\r\n";

  for (word = strtok(sentence, wordSeparators);
       i < max && word != 0;
       i++, word = strtok(0, wordSeparators)) {
    words[i] = word;
  }

  if (i < max) {
    words[i] = 0;
  }

  return i;
}

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

uint16_t big5touint(uint8_t h, uint8_t e)
{
	return ((((uint16_t)h << 8) ) | ((uint16_t)e));
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
	
	//output File
	FILE *fp;	
	fp = fopen(output_filename,"w+");
	map<uint16_t,vector<uint16_t>> zmap;
	File mapFile(map_filename, "r");
	char *line;
	Vocab all;
	while(line = mapFile.getline()){
		int line_len;
		char word[5];
		vector<uint16_t> vec;
		line_len = strlen(line);
		for(int i = 3;i<line_len;i+=3){
			memset(word, '\0',sizeof(word));
			strncpy(word,line+i*sizeof(char),2);
			vec.push_back(big5touint(word[0],word[1]));
		}
		memset(word, '\0',sizeof(word));
		strncpy(word,line,2);
		zmap[big5touint(word[0],word[1])]=vec;
	}
	mapFile.close();
	/*load LM*/
	Vocab voc;
	Ngram lm(voc, ngram_order);
	File lmFile(lm_filename, "r" );
        lm.read(lmFile);
        lmFile.close();
	File segFile(seg_filename, "r");
	while(line = segFile.getline()){
		//SRILM disambig output <s>X X X</
		char *words[500];
		char word[5];
		int count = parseWords(line,&(words[1]),500);
		words[0] = "<s>";
		words[count+1] = "</s>";
		count += 2;
		//build variable for viterbi 
		LogP deltap[200][1024] = {{0.0}};
		VocabIndex VidxGraph[200][1024];
		vector<uint16_t>::iterator it_i;
		vector<uint16_t> vec;
		int k,testsian[count];
		VocabIndex StrGraph[200][1024];
		int Backtrack[200][1024];
		int CandiNum[200];
		Prob p;
		VocabIndex empty[] = {Vocab_None};
		VocabIndex bi[] = {Vocab_None, Vocab_None};
		//initialize viterbi
		int size = 0;
		VocabIndex wid = voc.getIndex("<s>");
		LogP temp = lm.wordProb(wid, empty);//unigram
		if(temp== LogP_Zero){
			deltap[0][size] = -100;
		}else{
			deltap[0][size] = temp;
		}
		VidxGraph[0][size] = wid;
		StrGraph[0][size] = 0 ;
		Backtrack[0][size] = 0; 
		size++;
		CandiNum[0] = 1;
		//recursive solve viterbi
		for (int i = 1; i < count-1; i++) {
			size = 0;
			memset(word, '\0',sizeof(word));
			strncpy(word,words[i],2);
			vec = zmap.find(big5touint(word[0],word[1]))->second; 
			for(int why = 0;why<vec.size();why++){
				/*uint8_t vv[3] = {0};
				vv[0] = (uint8_t) (vec[why] >> 8);
				vv[1] = (uint8_t) (vec[why] & 0xFF);*/
				StrGraph[i][why] = vec[why];
			}
			testsian[i]= vec.size();
			for( k =0; k<vec.size();k++){
				uint8_t big5[3] ={0};
				big5[0] = (uint8_t) (vec[k] >> 8);
				big5[1] = (uint8_t) (vec[k] & 0xFF);
				VocabIndex wid = voc.getIndex((char*)big5);
				if(wid == Vocab_None){
					wid = voc.getIndex(Vocab_Unknown);
				}
				LogP maxp = LogP_Zero;
				for (int j = 0; j < CandiNum[i-1]; j++) {
						VocabIndex last;
						last = VidxGraph[i-1][j];
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
				size++;
			}
			
			CandiNum[i] = size;
		}
		testsian[0]= 1;
		testsian[count-1]= 1;
		size = 0;
		memset(word, '\0',sizeof(word));
		strncpy(word,words[count-1],4);
		wid = voc.getIndex(word);
		if(wid == Vocab_None){
			wid = voc.getIndex(Vocab_Unknown);
		}
		LogP maxp = LogP_Zero;
		for (int j = 0; j < CandiNum[count-2]; j++) {
				VocabIndex last;
				last = VidxGraph[count-2][j];
				if(last == Vocab_None){
					last = voc.getIndex(Vocab_Unknown);
				}						
				bi[0] = last;
				LogP logp = lm.wordProb(wid, bi); 
				if (logp == LogP_Zero ){
					logp = -100;
				} 
				logp += deltap[count-2][j];
				if (logp > maxp) {
					maxp = logp;
					Backtrack[count-1][size] = j;
				}
		}
		deltap[count-1][size] = maxp;
		VidxGraph[count-1][size] = wid;
		StrGraph[count-1][size] = 0;
		size++;
		
		CandiNum[count-1] = size;

		int finalmax = -1,j;
		LogP finalmaxp = LogP_Zero;
		for (j = 0; j < CandiNum[count-1]; j++) {
			if (deltap[count-1][j] > finalmaxp) {
				finalmaxp = deltap[count-1][j];
				finalmax = j;
			}
		}
		char AnsPath[maxWordLength][5];
		for(int i = 0;i<count;i++){
			for(int j = 0;j<testsian[i];j++){
				//printf("%u", StrGraph[i][j]);	
				uint8_t big5[3]={0};
				big5[0] = (uint8_t) (StrGraph[i][j] >> 8 );
				big5[1] = (uint8_t) (StrGraph[i][j] & 0xFF );
				//big5toutf8((char*)big5);
			}	
			//std::cout<<"--------------------------"<<endl;
		}
		for (int i = count-1; i > 0; i--) {
			//uint8_t big5[3]={0};
			AnsPath[i][0] = (char) (StrGraph[i][finalmax] >> 8 );
			AnsPath[i][1] = (char) (StrGraph[i][finalmax] & 0xFF );
			AnsPath[i][2] = (char) 0; 
			//printf("%d ",finalmax);
			//big5toutf8((char*)big5);
			//AnsPath[i][] = (char*)big5;
			//AnsPath[i] = StrGraph[i][finalmax];
			//AnsPath[i] = voc.getWord(VidxGraph[i][finalmax]);
			finalmax = Backtrack[i][finalmax];
		}
		AnsPath[0][0] = '<';
		AnsPath[0][1] = 's';
		AnsPath[0][2] = '>';
		AnsPath[0][3] = '\0';
		AnsPath[count-1][0] = '<';
		AnsPath[count-1][1] = '/';
		AnsPath[count-1][2] = 's';
		AnsPath[count-1][3] = '>';
		AnsPath[count-1][4] = '\0';
		// Print the Answer Path
		for (int i = 0; i < count; i++){
			fprintf(fp,"%s%s", AnsPath[i], (i == count-1)? "\n": " ");
		}
	}
	fclose(fp);
	return 0;
}
