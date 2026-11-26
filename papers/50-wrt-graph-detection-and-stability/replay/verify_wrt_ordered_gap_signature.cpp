#include<bits/stdc++.h>
using namespace std;
using ll=long long;
int M, na, nb;
struct Key {int a,k; string sign; bool operator<(Key const&o)const {return tie(a,k,sign)<tie(o.a,o.k,o.sign);}};
map<Key,ll> diff;
int sg,a,cap;
vector<pair<int,int>> es;
int positions[32], fmask[30], fcount;
vector<int> chosen;
void record(int unionmask) {
 int k=__builtin_popcount((unsigned)unionmask);
 int v[3]={},z=0;for(int m=unionmask;m;m&=m-1)v[z++]=m&-m;
 int p[3]={0,1,2};string best;
 do {
   string s;
   for(int c:chosen) {
     int mask=abs(c), w=0; for(int j=0;j<k;j++)if(v[p[j]]&mask)w|=1<<j;
     s.push_back((c>0?'a':'A')+w-1);
   }
   sort(s.begin(),s.end());
   if(best.empty() || s<best)best=s;
 }while(next_permutation(p,p+k));
 diff[{a,k,best}]+=sg;
}
void dfs(int last,int mask) {
 record(mask);
 for(int j=last;j<fcount;j++)if(__builtin_popcount((unsigned)(mask|fmask[j]))<=cap) {
    chosen.push_back(fmask[j]);dfs(j+1,mask|fmask[j]);chosen.pop_back();
 }
}
void getorder(int n,bool bpin[]) {
 int min=0;fcount=0;chosen.clear();
 for (auto [i,j]:es) {
    if(bpin[i] && bpin[j])continue;
    int p=positions[i],q=positions[j];int m=((1<<max(p,q))-1)^((1<<std::min(p,q))-1);
    if(p>q) {
      min|=m;
      if(__builtin_popcount((unsigned)min)>cap)return;
      chosen.push_back(-m);
    } else if(__builtin_popcount((unsigned)m)<=cap) fmask[fcount++]=m;
 }
 dfs(0,min);
}
void graph(vector<int> r,int sign) {
 sg=sign;es.clear();
 for(int i=0;i<na;i++)for(int j=0;j<nb;j++)if(r[i]>>j&1)es.push_back({i,na+j});
 for(int b=0;b<(1<<M);b++) {
   a=__builtin_popcount((unsigned)b);if(a>3)continue;cap=3-a;
   vector<int> pp;bool bpin[M];int n=M-a;
   for(int i=0;i<M;i++){bpin[i]=b>>i&1;if(!bpin[i])pp.push_back(i);else positions[i]=i<na?0:n+1;}
   do {
      for(int i=0;i<n;i++)positions[pp[i]]=i+1;
      getorder(n,bpin);
   }while(next_permutation(pp.begin(),pp.end()));
 }
}
int main(int argc, char** argv){
 na=stoi(argv[1]);nb=stoi(argv[2]);M=na+nb; vector<int> p,q;for(int i=0;i<na;i++)p.push_back(stoi(argv[i+3])); graph(p,1);if(argc>3+na){for(int i=0;i<na;i++)q.push_back(stoi(argv[i+3+na]));graph(q,-1);}
 for(auto [k,c]:diff)if(c)cout<<k.a<<' '<<k.k<<' '<<(k.sign.empty()? "-":k.sign)<<' '<<c<<'\n';
}
