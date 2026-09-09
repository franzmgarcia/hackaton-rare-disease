// Bounded-memory targeted FASTQ pair extraction. Retains both mates if either
// contains a reference/alternate 31-mer seed. This is candidate retrieval,
// NOT alignment, phase calling, or an exhaustive error-tolerant mapper.
#include <zlib.h>
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_set>
#include <vector>
using namespace std;
const uint64_t MASK=(1ULL<<62)-1;
int base(char c){switch(c){case 'A':case 'a':return 0;case 'C':case 'c':return 1;case 'G':case 'g':return 2;case 'T':case 't':return 3;default:return -1;}}
uint64_t code(const string&s){uint64_t f=0,r=0;for(char c:s){int b=base(c);if(b<0)throw runtime_error("invalid seed");f=(f<<2)|b;r=(r>>2)|((uint64_t)(3-b)<<60);}return min(f,r);}
bool hit(const string&s,const unordered_set<uint64_t>&seeds){uint64_t f=0,r=0;int valid=0;for(char c:s){int b=base(c);if(b<0){valid=0;f=r=0;continue;}f=((f<<2)|b)&MASK;r=(r>>2)|((uint64_t)(3-b)<<60);if(++valid>=31 && seeds.count(min(f,r)))return true;}return false;}
bool line(gzFile f,string&s){s.clear();char b[4096];while(gzgets(f,b,sizeof(b))){s+=b;if(s.back()=='\n'){s.pop_back();return true;}}if(!s.empty())throw runtime_error("truncated FASTQ line");int e;gzerror(f,&e);if(e!=Z_OK&&e!=Z_STREAM_END)throw runtime_error("gzip failure");return false;}
bool record(gzFile f,vector<string>&r){if(!line(f,r[0]))return false;for(int j=1;j<4;j++)if(!line(f,r[j]))throw runtime_error("truncated FASTQ");if(r[0][0]!='@'||r[2][0]!='+'||r[1].size()!=r[3].size())throw runtime_error("invalid FASTQ");return true;}
string id(string s){s=s.substr(0,s.find(' '));if(s.size()>2&&s[s.size()-2]=='/'&&(s.back()=='1'||s.back()=='2'))s.resize(s.size()-2);return s;}
int main(int argc,char**argv){try{if(argc!=6)throw runtime_error("usage seeds R1.gz R2.gz out1.fq out2.fq");unordered_set<uint64_t>seeds;ifstream sf(argv[1]);string s;while(sf>>s){if(s.size()!=31)throw runtime_error("seed length");seeds.insert(code(s));}if(seeds.empty())throw runtime_error("no seeds");gzFile a=gzopen(argv[2],"rb"),b=gzopen(argv[3],"rb");if(!a||!b)throw runtime_error("input open");gzbuffer(a,1<<20);gzbuffer(b,1<<20);ofstream oa(argv[4]),ob(argv[5]);if(!oa||!ob)throw runtime_error("output open");vector<string>x(4),y(4);uint64_t total=0,kept=0,minlen=100000,maxlen=0;while(record(a,x)){if(!record(b,y)||id(x[0])!=id(y[0]))throw runtime_error("mates out of sync");total++;for(auto*r:{&x,&y}){minlen=min(minlen,(uint64_t)(*r)[1].size());maxlen=max(maxlen,(uint64_t)(*r)[1].size());}if(hit(x[1],seeds)||hit(y[1],seeds)){for(int j=0;j<4;j++){oa<<x[j]<<'\n';ob<<y[j]<<'\n';}kept++;}if(total%10000000==0)cerr<<"pairs="<<total<<" retained="<<kept<<endl;}if(record(b,y))throw runtime_error("R2 extra records");gzclose(a);gzclose(b);cout<<"{\"pairs_scanned\":"<<total<<",\"pairs_retained\":"<<kept<<",\"read_length_min\":"<<minlen<<",\"read_length_max\":"<<maxlen<<",\"seeds\":"<<seeds.size()<<"}"<<endl;return 0;}catch(const exception&e){cerr<<e.what()<<endl;return 1;}}
