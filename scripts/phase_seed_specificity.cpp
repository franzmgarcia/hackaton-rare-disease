// Count seed occurrences across the full reference, saturating at two.
#define main extraction_main
#include "phase_extract.cpp"
#undef main
#include <unordered_map>
int main(int argc,char**argv){try{
 if(argc!=5)throw runtime_error("usage seeds genome.fa.gz selected_seeds counts");
 ifstream in(argv[1]);string s;vector<string>strings;unordered_map<uint64_t,int>counts;
 while(in>>s){strings.push_back(s);counts[code(s)]=0;}
 gzFile f=gzopen(argv[2],"rb");if(!f)throw runtime_error("open genome");gzbuffer(f,1<<20);
 uint64_t forward=0,reverse=0,total=0;int valid=0;
 while(line(f,s)){if(s[0]=='>'){valid=0;forward=reverse=0;continue;}
  for(char c:s){int b=base(c);total++;if(b<0){valid=0;forward=reverse=0;continue;}
   forward=((forward<<2)|b)&MASK;reverse=(reverse>>2)|((uint64_t)(3-b)<<60);
   if(++valid>=31){auto it=counts.find(min(forward,reverse));if(it!=counts.end()&&it->second<2)it->second++;}
  }
 }
 gzclose(f);ofstream selected(argv[3]),audit(argv[4]);int n=0;
 for(auto&seed:strings){int c=counts[code(seed)];audit<<seed<<'\t'<<c<<'\n';if(c<=1){selected<<seed<<'\n';n++;}}
 cout<<"{\"genome_bases_scanned\":"<<total<<",\"input_seeds\":"<<strings.size()<<",\"selected_seeds\":"<<n<<"}"<<endl;
 return 0;
 }catch(const exception&e){cerr<<e.what()<<endl;return 1;}}
